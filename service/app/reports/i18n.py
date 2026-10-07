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
