# Docker イメージのスリム化・ビルド高速化 設計書

- 日付: 2026-06-03
- 対象: `dfast_vrl` の Docker / Singularity ビルド
- ステータス: 設計承認済み（実装計画へ）

## 目的

`nigyta/dfast_vrl` コンテナの **イメージサイズ削減** と **ビルド時間短縮** を、
挙動を変えずに両立する。スリム化と高速化は両方バランス良く狙う。

## 背景・現状の問題点

現行 `Dockerfile`（`ubuntu` ベース）の主なコスト要因:

**ビルド時間**
- `vadr-install.sh` が VADR を**ソースからコンパイル**（infernal / Bio-Easel /
  HMMER / BLAST 等）。単一で最大のコスト。
- `mamba install`（`snpeff=5.0` が JRE を引き込む + 依存解決）も重い。

**イメージサイズ**
1. `build-essential` / `git` / `autoconf` がランタイムにも残る（ビルド時のみ必要）。
2. `mamba clean` 未実施 → conda パッケージキャッシュが数百 MB 残存。
3. apt の `rm -rf /var/lib/apt/lists/*` が実コードに無い（コメントのみ）。
4. snpEff 同梱 JRE が大きい。
5. `.dockerignore` が無い → `OUT_*` / `examples` / `.git` などが毎回ビルド
   コンテキストとして送られる（CI は `context: .`）。

## コード調査で確定した事実

| 項目 | 結論 |
|---|---|
| `pandas` | コード内に参照ゼロ → 明示インストールから**削除可**（dr_tools が必要とすれば pip が依存として自動導入） |
| `biopython` | 10 ファイルで使用 → 必須（`pip` 経由） |
| `mafft` | `dfv/detect_variants.py` で `mafft` コマンド呼び出し → 必須 |
| `snpEff` | `dfv/detect_variants.py:140` で `snpEff -c ...` をコマンド実行。設定 `refs/snpeff/snpEff.config`（`data.dir=./data/`）、DB `nigvrl`/`nigvrl2` は `refs/snpeff/data/` にコミット済 |
| `.bin` 互換 | `snpEffectPredictor.bin` は snpEff **5.0** でビルド済 → snpEff は 5.0 にピン必須 |
| `dr_tools` | `dfv/common.py:mss2json` が `from dr_tools import drt_ann2json`（try/except で握り潰し）→ `pip`(git) 必要 |
| VADR models | コンテナに含めず、ホストの `VADRMODELDIR=/vadr_models` をマウント（現状方針どおり） |
| shebang | 3 スクリプトとも `#!/usr/bin/env python` → `python` が PATH に必要 |

### ベースイメージ調査（staphb/vadr）

- **どの staphb/vadr イメージにも Java は含まれない**（openjdk/jre のインストール無し）
  → snpEff 用に **JRE の追加が必須**。
- staphb/vadr は `VADRMODELDIR=/opt/vadr/vadr-models` を設定（本プロジェクトは
  `/vadr_models` にホストマウント）→ base の ENV を継承しつつ `VADRMODELDIR` のみ上書き。
- registry 圧縮サイズ: `1.6.4`=542MB / **`1.7-slim`=583MB** / `1.7`=840MB。
- 採用タグ: **`staphb/vadr:1.7-slim`**（ubuntu:noble ベースの軽量版）。

### snpEff 5.0 入手

- `https://snpeff-public.s3.amazonaws.com/versions/snpEff_v5_0_core.zip`
  → HTTP 200 / 45.8MB を確認済み。
- Java: snpEff 5.0 は Java 11〜17 で安定（21 では旧 5.0 で不具合報告あり）。
  → `openjdk-17-jre-headless` を採用。検証で不可なら 11 系へフォールバック。

## 決定事項（ユーザー承認済み）

- VADR 導入: **prebuilt `staphb/vadr:1.7-slim` を base に**（コンパイル省略）。
- VADR モデルは**コンテナに含めず**ホストマウント（`/vadr_models`）。
- **conda/mamba を全廃**し、各依存を最小手段で導入。
- 優先度: スリム化とビルド時間短縮を**両方バランス良く**。
- **既存 `Dockerfile` / `Singularity.def` は変更せず温存**し、`*.slim` を新規作成。
- Docker と Singularity を**同方針で同期**。
- pip は `--break-system-packages`（PEP668 回避、単一用途コンテナのため許容）。

## 成果物（新規ファイル）

- `Dockerfile.slim`
- `Singularity.slim.def`（同方針ミラー）
- `.dockerignore`

既存の `Dockerfile` / `Singularity.def` は今回**変更しない**。

## 設計詳細

### ベース・環境変数

```dockerfile
FROM --platform=linux/amd64 staphb/vadr:1.7-slim
```

- base が設定する VADR 系 ENV（`VADRINSTALLDIR`, `VADRSCRIPTSDIR`,
  `VADRMINISCRIPTSDIR`, `PERL5LIB`, VADR の `PATH` 等）は**継承する**。
- 本プロジェクト固有の上書きのみ宣言:
  - `VADRMODELDIR=/vadr_models`（ホストマウント先）
  - `LC_ALL=C`
  - snpEff ラッパー設置先 `/usr/local/bin` は既定 PATH 内のため PATH 追記不要。

### 依存の置換

| 依存 | 置換方法 |
|---|---|
| mafft | `apt-get install -y --no-install-recommends mafft` |
| Java | `apt-get install -y --no-install-recommends openjdk-17-jre-headless` |
| snpEff 5.0 | S3 zip を `/opt/snpeff` に展開し、`/usr/local/bin/snpEff` ラッパー作成 |
| python | `apt-get install -y --no-install-recommends python3 python3-pip python-is-python3` |
| biopython | `pip install --no-cache-dir --break-system-packages biopython==1.84` |
| dr_tools | `pip install --no-cache-dir --break-system-packages "git+https://github.com/ddbj/dr_tools.git"` |
| pandas | **入れない** |

snpEff ラッパー（例）:

```sh
#!/bin/sh
exec java -jar /opt/snpeff/snpEff/snpEff.jar "$@"
```

（zip 展開後の実際のディレクトリ構成・jar パスは実装時に確認して合わせる。
`snpEff -c <config>` で起動するため、ラッパーが引数を透過すれば既存コードは無改修。）

### スリム化施策

1. `.dockerignore` を新規作成し、ビルドに不要なものを除外:
   `OUT*`, `examples`, `.git`, `dev`, `tests`, `cox1_test`, `*_test`,
   `.pytest_cache`, `.devcontainer`, `.DS_Store`, `docs` など。
   （本体は `git clone` で取得するため、コンテキストは小さくてよい）
2. base 同梱の不要モデルが残っていれば `rm -rf /opt/vadr/vadr-models*`
   （ホストマウント方針のため不要）。
3. 各 `RUN` 末で `apt-get clean && rm -rf /var/lib/apt/lists/*`、
   `pip --no-cache-dir`、snpEff zip の削除。
4. RUN レイヤーを意味単位で集約。

### dfast_vrl 本体導入（現行踏襲）

```dockerfile
ARG INCREMENT_THIS_TO_DISABLE_CACHE_BELOW_THIS_LINE=1
RUN git clone https://github.com/nigyta/dfast_vrl.git && \
    ln -s /dfast_vrl/dfast_vrl /usr/bin/dfast_vrl && \
    ln -s /dfast_vrl/vadr2mss.py /usr/bin/vadr2mss.py
WORKDIR /data
RUN dfast_vrl --version
```

### Singularity.slim.def

`Bootstrap: docker` / `From: staphb/vadr:1.7-slim` とし、`%post` を
Dockerfile.slim と同じ手順（apt 依存、snpEff、pip、本体 clone）で構成。
`%environment` は VADR base の前提に合わせ `VADRMODELDIR=/vadr_models` 等を export。

## リスクと検証

### リスク

- **VADR 1.6.4 → 1.7 のモデル互換性**: 現ピン（SCOV2 1.3-2 / MPXV 1.4.2-1 /
  RSV 1.5-2 / COX1 1.2-1）が VADR 1.7 で動くか未確認。動かない場合はモデル版の
  更新が必要。→ 検証で確認。
- **snpEff 5.0 と Java 17 の組合せ**: `.bin` 読み込み・実行が通るか。
  不可なら JRE を 11 系へ変更。→ 検証で確認。

### 検証（受け入れ条件）

1. `docker build -f Dockerfile.slim -t dfast_vrl:slim .` が成功する。
2. `docker run --rm dfast_vrl:slim dfast_vrl --version` が動く。
3. SARS-CoV-2 パイプライン 1 本（`dfast_vrl`、mafft+snpEff 経路）を
   `examples/` の入力で実行し完走する（モデルはマウント）。
4. `docker images` で**旧 `Dockerfile` ビルドとサイズを数値比較**し削減を確認。
5. 可能なら `vadr2mss.py`（汎用 VADR 経路）も 1 本実行（モデルがあれば）。

## スコープ外（今回触れない）

- CI（`.github/workflows/update_container.yaml`）の Dockerfile 切替。
  slim を用意したのち別途相談する。
- 既存 `Dockerfile` / `Singularity.def` の変更。
- CLI スクリプトのバージョン文字列等、Docker と無関係な変更。
