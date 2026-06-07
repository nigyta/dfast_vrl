# staphb/vadr:1.7-slim — フォールバック用ビルド資料

`Dockerfile` は `staphb/vadr:1.7-slim` の元 Dockerfile（StaPH-B/docker-builds）の
バイト等価コピー。本プロジェクトの `Dockerfile` / `Singularity.def` は
`staphb/vadr:1.7-slim` を base にしているため、**万一この base イメージを取得できない
場合に、ここから base を自前ビルドできる**ように保存している。

- 取得元: <https://github.com/StaPH-B/docker-builds/blob/master/build-files/vadr/1.7-slim/Dockerfile>
- raw: <https://raw.githubusercontent.com/StaPH-B/docker-builds/master/build-files/vadr/1.7-slim/Dockerfile>
- 取得日: 2026-06-08
- 公開イメージ `1.7-slim` は **`app` ステージ**に相当（`test` ステージは
  sarscov2 モデルDL＋動作テスト用で、公開イメージには含まれない）。

## 自前ビルド

```bash
# このディレクトリで（app ステージのみ。test ステージはビルドしない）
docker build --platform=linux/amd64 --target app -t staphb/vadr:1.7-slim .
```

- **ネットワーク必須**: VADR ソースを GitHub releases から取得し、`cpan` で
  `Mozilla::CA` を導入する。
- **ローカルのビルドコンテキストは不要**: `app` ステージは `COPY` を使わない
  （`docker build ... - < Dockerfile` のように stdin からでもビルド可能）。
- **注意**: この Dockerfile は `vadr-install.sh` の特定行を `sed -i '228,278d'` で
  削除してデフォルトモデルのダウンロードを省く実装。VADR タグ側の更新で行番号が
  ずれると壊れうるため、これは「取得時点（2026-06-08）の記録」として保存している。

ビルド後は本プロジェクトの `Dockerfile` がそのまま
`FROM --platform=linux/amd64 staphb/vadr:1.7-slim` で利用できる。
