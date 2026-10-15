"""Report i18n string catalogs (EN/FR/HT).

Source of truth for every human-readable string in the rendered report.
English is the reference language. French and Haitian Creole strings are
DRAFT — they carry a TRANSLATOR-REVIEW marker until a native speaker signs
off (see render: non-English strings are wrapped in an HTML comment AND the
report shows a visible review banner for fr/ht).

Canonical strings (verbatim, from sentinel-liability-guardrails.md):
- RUO_LABEL_* : "For Research Use Only. Not for use in diagnostic procedures."
- DISCLAIMER_* : the extended surveillance disclaimer (all three languages).
These are the ONLY strings permitted to contain otherwise-forbidden words
(diagnose/treat/patient...) — the banned-phrase CI test strips them first.
"""

from __future__ import annotations

# Marker wrapped around every non-English string at render time.
TRANSLATOR_REVIEW = "<!-- TRANSLATOR-REVIEW -->"

RUO_LABEL = {
    "en": "For Research Use Only. Not for use in diagnostic procedures.",
    # DRAFT — native-speaker review required.
    "fr": "Pour usage en recherche uniquement. Ne pas utiliser dans des procédures diagnostiques.",
    # DRAFT — native-speaker review required.
    "ht": "Pou itilizasyon rechèch sèlman. Pa itilize nan pwosedi dyagnostik.",
}

DISCLAIMER = {
    "en": (
        "This report was generated for public-health surveillance and research "
        "purposes only. Results are computational predictions derived from "
        "genomic sequence data and must be confirmed by phenotypic laboratory "
        "testing. Do not use this report to diagnose or treat any individual "
        "patient. Antimicrobial resistance markers reported here are genetic "
        "observations, not susceptibility results."
    ),
    # DRAFT — from the guardrails doc; native-speaker review required.
    "fr": (
        "Ce rapport est établi uniquement à des fins de surveillance de santé "
        "publique et de recherche. Les résultats sont des prédictions "
        "informatiques issues de données de séquençage génomique et doivent "
        "être confirmés par des tests phénotypiques en laboratoire. Ne pas "
        "utiliser ce rapport pour diagnostiquer ou traiter un patient. Les "
        "marqueurs de résistance aux antimicrobiens rapportés ici sont des "
        "observations génétiques, pas des résultats de sensibilité."
    ),
    # DRAFT — from the guardrails doc; native-speaker review required.
    "ht": (
        "Rapò sa a fèt pou siveyans sante piblik ak rechèch sèlman. Rezilta yo "
        "se prediksyon òdinatè ki soti nan done sekans jenomik; fòk yo konfime "
        "yo ak tès laboratwa (fenotipik). Pa itilize rapò sa a pou dyagnostike "
        "oswa trete okenn pasyan. Makè rezistans antimikwòb yo rapòte isit la "
        "se obsèvasyon jenetik, se pa rezilta sansiblite."
    ),
}

STRINGS: dict[str, dict[str, str]] = {
    "en": {
        "title": "Genomic surveillance report",
        "ruo_banner": "For Research Use Only. Not for use in diagnostic procedures.",
        "qc_verdict": "QC verdict",
        "qc_pass": "PASS",
        "qc_fail": "FAIL",
        "not_interpretable": "Results below are not interpretable.",
        "qc_reasons": "QC reasons",
        "section_sample": "Sample",
        "section_species": "Species identification",
        "section_lineage": "Lineage",
        "section_loci": "Surveillance loci",
        "section_amr": "AMR marker detection",
        "section_provenance": "Provenance",
        "label_job": "Job / sample ID",
        "label_org": "Organization",
        "label_uploaded": "Uploaded",
        "label_completed": "Completed",
        "label_pipeline": "Pipeline version",
        "label_reference": "Reference",
        "label_platform": "Platform",
        "label_basecaller": "Basecaller model",
        "label_tool_versions": "Tool and database versions",
        "label_species": "Species",
        "label_species_reads": "V. cholerae reads",
        "label_species_fraction": "Fraction of classified reads",
        "label_snps": "SNP distance vs 7PET reference",
        "label_mlst": "MLST sequence type",
        "label_vibecheck": "Vibecheck lineage",
        "label_mean_depth": "Mean depth",
        "label_breadth": "Genome breadth covered",
        "locus": "Locus",
        "call": "Call",
        "present": "present",
        "partial": "partial",
        "absent": "absent",
        "gene": "Gene",
        "detected": "detected",
        "not_detected": "not detected",
        "amr_disclaimer": (
            "Gene detected / not detected only. These are genetic observations, "
            "not susceptibility results."
        ),
        "footer_note": (
            "Interpretation aid for public-health surveillance and research. "
            "Not a medical device."
        ),
        "generated_at": "Report generated",
        "no_data": "—",
        "section_confidence": "Confidence",
        "confidence_high": "High",
        "confidence_provisional": "Provisional",
        "confidence_low": "Low — do not act",
        "confidence_operational_note": (
            "Operational confidence tier, not a published standard."
        ),
        "confidence_reasons": "Reasons",
        "reason_qc_fail": "QC failed — results are not interpretable.",
        "reason_depth": (
            "mean depth {depth}x is below the 30x high-confidence threshold"
        ),
        "reason_called": (
            "{called}% of the consensus is called, "
            "below the 95% high-confidence threshold"
        ),
        "reason_contamination": (
            "{pct}% of reads are non-V. cholerae, "
            "above the 1% high-confidence threshold"
        ),
        "section_limitations": "Limitations",
        "limitation_1": (
            "SNP distances reflect recent common ancestry, not who infected "
            "whom — field epidemiology is still needed to infer transmission."
        ),
        "limitation_2": (
            "No fixed SNP cutoff defines an outbreak. The ≤5 SNP / 14-day "
            "alert rule is an operational choice, not a biological threshold."
        ),
        "limitation_3": (
            "Recombination is not masked in this screen. In Haiti, small SNP "
            "differences can reflect environmental adaptation or bottlenecks "
            "rather than transmission."
        ),
        "limitation_4": (
            "Mapping is reference-based: insertions and genes absent from the "
            "2010EL-1786 reference are invisible here. Mobile elements and "
            "AMR outside the core genome need assembly-based checks."
        ),
        "limitation_5": (
            "Gene presence is not phenotype. AMR calls below are gene "
            "detected / not detected only."
        ),
    },
    "fr": {
        # DRAFT — every string below requires native-speaker review.
        "title": "Rapport de surveillance génomique",
        "ruo_banner": "Pour usage en recherche uniquement. Ne pas utiliser dans des procédures diagnostiques.",
        "qc_verdict": "Verdict du contrôle qualité",
        "qc_pass": "RÉUSSI",
        "qc_fail": "ÉCHEC",
        "not_interpretable": "Les résultats ci-dessous ne sont pas interprétables.",
        "qc_reasons": "Motifs du contrôle qualité",
        "section_sample": "Échantillon",
        "section_species": "Identification de l'espèce",
        "section_lineage": "Lignage",
        "section_loci": "Loci de surveillance",
        "section_amr": "Détection de marqueurs de résistance",
        "section_provenance": "Provenance",
        "label_job": "ID du job / échantillon",
        "label_org": "Organisation",
        "label_uploaded": "Téléversé le",
        "label_completed": "Terminé le",
        "label_pipeline": "Version du pipeline",
        "label_reference": "Référence",
        "label_platform": "Plateforme",
        "label_basecaller": "Modèle de basecalling",
        "label_tool_versions": "Versions des outils et bases de données",
        "label_species": "Espèce",
        "label_species_reads": "Lectures V. cholerae",
        "label_species_fraction": "Fraction des lectures classées",
        "label_snps": "Distance SNP vs référence 7PET",
        "label_mlst": "Type de séquence MLST",
        "label_vibecheck": "Lignage Vibecheck",
        "label_mean_depth": "Profondeur moyenne",
        "label_breadth": "Couverture du génome",
        "locus": "Locus",
        "call": "Appel",
        "present": "présent",
        "partial": "partiel",
        "absent": "absent",
        "gene": "Gène",
        "detected": "détecté",
        "not_detected": "non détecté",
        "amr_disclaimer": (
            "Gène détecté / non détecté uniquement. Il s'agit d'observations "
            "génétiques, pas de résultats de sensibilité."
        ),
        "footer_note": (
            "Outil d'aide à l'interprétation pour la surveillance de santé "
            "publique et la recherche. Pas un dispositif médical."
        ),
        "generated_at": "Rapport généré le",
        "no_data": "—",
        "section_confidence": "Confiance",
        "confidence_high": "Élevé",
        "confidence_provisional": "Provisoire",
        "confidence_low": "Faible — ne pas agir",
        "confidence_operational_note": (
            "Niveau de confiance opérationnel, pas une norme publiée."
        ),
        "confidence_reasons": "Motifs",
        "reason_qc_fail": (
            "le contrôle qualité a échoué — résultats non interprétables."
        ),
        "reason_depth": (
            "la profondeur moyenne de {depth}x est sous le seuil de 30x"
        ),
        "reason_called": (
            "{called} % du consensus est appelé, sous le seuil de 95 %"
        ),
        "reason_contamination": (
            "{pct} % des lectures ne sont pas V. cholerae, "
            "au-dessus du seuil de 1 %"
        ),
        "section_limitations": "Limites",
        "limitation_1": (
            "Les distances en SNP reflètent une ascendance commune récente, "
            "pas l'identité de l'infecteur — l'épidémiologie de terrain "
            "reste nécessaire pour inférer la transmission."
        ),
        "limitation_2": (
            "Aucun seuil fixe de SNP ne définit une flambée. La règle "
            "d'alerte ≤5 SNP / 14 jours est un choix opérationnel, pas un "
            "seuil biologique."
        ),
        "limitation_3": (
            "La recombinaison n'est pas masquée dans cet écran. En Haïti, "
            "de petites différences de SNP peuvent refléter une adaptation "
            "environnementale ou des goulots d'étranglement plutôt qu'une "
            "transmission."
        ),
        "limitation_4": (
            "L'alignement est basé sur une référence : les insertions et "
            "les gènes absents de la référence 2010EL-1786 sont invisibles "
            "ici. Les éléments mobiles et la résistance aux antimicrobiens "
            "hors du génome cœur nécessitent des vérifications par assemblage."
        ),
        "limitation_5": (
            "La présence d'un gène ne prédit pas le phénotype. Les appels "
            "ci-dessous signifient « gène détecté / non détecté » uniquement."
        ),
    },
    "ht": {
        # DRAFT — every string below requires native-speaker review.
        "title": "Rapò siveyans jenomik",
        "ruo_banner": "Pou itilizasyon rechèch sèlman. Pa itilize nan pwosedi dyagnostik.",
        "qc_verdict": "Verdik kontwòl kalite",
        "qc_pass": "PASE",
        "qc_fail": "ECHWE",
        "not_interpretable": "Rezilta ki anba yo pa ka entèprete.",
        "qc_reasons": "Rezon kontwòl kalite",
        "section_sample": "Echantiyon",
        "section_species": "Idantifikasyon espès",
        "section_lineage": "Liyaj",
        "section_loci": "Loci siveyans",
        "section_amr": "Deteksyon makè rezistans",
        "section_provenance": "Pwovnans",
        "label_job": "ID job / echantiyon",
        "label_org": "Òganizasyon",
        "label_uploaded": "Te voye",
        "label_completed": "Te fini",
        "label_pipeline": "Vèsyon pipeline",
        "label_reference": "Referans",
        "label_platform": "Platfòm",
        "label_basecaller": "Modèl basecalling",
        "label_tool_versions": "Vèsyon zouti ak baz done",
        "label_species": "Espès",
        "label_species_reads": "Lekti V. cholerae",
        "label_species_fraction": "Fraksyon lekti klase yo",
        "label_snps": "Distans SNP ak referans 7PET",
        "label_mlst": "Tip sekans MLST",
        "label_vibecheck": "Liyaj Vibecheck",
        "label_mean_depth": "Pwofondè mwayèn",
        "label_breadth": "Kouvèti jenom",
        "locus": "Locus",
        "call": "Apèl",
        "present": "prezan",
        "partial": "pasyèl",
        "absent": "absan",
        "gene": "Jèn",
        "detected": "detekte",
        "not_detected": "pa detekte",
        "amr_disclaimer": (
            "Jèn detekte / pa detekte sèlman. Se obsèvasyon jenetik, "
            "se pa rezilta sansiblite."
        ),
        "footer_note": (
            "Zouti èd pou entèprete pou siveyans sante piblik ak rechèch. "
            "Se pa yon aparèy medikal."
        ),
        "generated_at": "Rapò a te pwodui",
        "no_data": "—",
        "section_confidence": "Konfyans",
        "confidence_high": "Wo",
        "confidence_provisional": "Pwovizwa",
        "confidence_low": "Fèb — pa aji",
        "confidence_operational_note": (
            "Nivo konfyans operasyonèl, se pa yon estanda pibliye."
        ),
        "confidence_reasons": "Rezon",
        "reason_qc_fail": (
            "kontwòl kalite a echwe — rezilta yo pa entèprete."
        ),
        "reason_depth": (
            "pwofondè mwayèn {depth}x la pi ba pase papòt 30x la"
        ),
        "reason_called": (
            "{called} % konsensis la rele, pi ba pase papòt 95 % la"
        ),
        "reason_contamination": (
            "{pct} % lekti yo se pa V. cholerae, pi wo pase papòt 1 % la"
        ),
        "section_limitations": "Limit",
        "limitation_1": (
            "Distans SNP yo montre yon zansèt komen resan, pa ki moun ki "
            "bay ki moun maladi a — epidemyoloji teren toujou nesesè pou "
            "konprann transmisyon an."
        ),
        "limitation_2": (
            "Okenn chif SNP fiks pa defini yon epidemi. Règ alèt ≤5 SNP / "
            "14 jou a se yon chwa operasyonèl, se pa yon papòt byolojik."
        ),
        "limitation_3": (
            "Rekonbinasyon pa kache nan egzamen sa a. An Ayiti, ti diferans "
            "SNP yo ka reflete adaptasyon anviwònman oswa peryòd blokaj, "
            "pa transmisyon."
        ),
        "limitation_4": (
            "Aliyman an fèt pa referans: ensèsyon ak jèn ki pa nan referans "
            "2010EL-1786 an pa parèt isit la. Eleman mobil ak rezistans "
            "antimikwòb deyò jenom nwayo a bezwen verifikasyon pa asanblaj."
        ),
        "limitation_5": (
            "Prezans yon jèn pa vle di fenotip. Apèl yo pi ba a vle di "
            "« jèn detekte / pa detekte » sèlman."
        ),
    },
}

# Visible banner shown at the top of every non-English report.
REVIEW_BANNERS = {
    "fr": (
        "Traduction en cours : les libellés en français sont en attente de "
        "révision par un locuteur natif."
    ),
    "ht": (
        "Tradiksyon an pwovizwa: tout tèks kreyòl yo ap tann revizyon "
        "yon moun ki pale kreyòl kòm lang matènèl."
    ),
}


def get_strings(lang: str) -> tuple[dict[str, str], bool]:
    """(strings, needs_review). Unknown langs fall back to English."""
    if lang in STRINGS:
        return STRINGS[lang], lang != "en"
    return STRINGS["en"], False


def mark_review(text: str, needs_review: bool) -> str:
    """Wrap a non-English string in the TRANSLATOR-REVIEW marker."""
    if needs_review:
        return f"{TRANSLATOR_REVIEW}{text}"
    return text
