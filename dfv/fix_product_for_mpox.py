
dict_product = {
    "A32.5L": "A32.5L protein",
    "M137R": "M137R protein",
    "MPXVgp154": "hypothetical protein",
    "A5L protein-like": "A5L-like protein",
    "Brix domain protein": "Brix domain-containing protein",
    "DUF1029 domain protein": "DUF1029 domain-containing protein",
    "EGF-like domain protein": "EGF-like domain-containing protein",
    "NUDIX domain protein": "NUDIX domain-containing protein",
    "Profilin domain protein": "Profilin domain-containing protein",
    "RNA polymerase-associate": "RNA polymerase-associate protein",
    "Crescent membrane and immature virion formatio": "Crescent membrane and immature virion formation protein",
    "MPXVgp007": "D4L protein",
    "MPXVgp021": "D18L protein",
    "MPV-Z-N2R": "MPV-Z-N2R protein",
    "A-type inclusion protein (Cop-A25L); A26L": "A-type inclusion protein",
    "A-type inclusion protein (Cop-A25L); A27L": "A-type inclusion protein",
    "BTB Kelch-domain containing protein CRL complex (Cop-A55R); B1R kelch-like": "BTB Kelch-domain containing protein CRL complex",
    "Host range protein; D2L": "Host range protein",
    "IL-1 beta inhibitor (Cop-B17L); B15L": "IL-1 beta inhibitor",
    "NMDA receptor-like protein; R1R": "NMDA receptor-like protein",
    "kelch-like protein (EV-M-167); B18R kelch-like fragment": "kelch-like protein"
}

def fix_product_for_mpox(original_product):
    if original_product in dict_product:
        return dict_product[original_product]
    else:
        return original_product
