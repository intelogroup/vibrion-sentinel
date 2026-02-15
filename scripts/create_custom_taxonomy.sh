#!/bin/bash
# Create custom taxonomy for Kraken2 Haiti database
# Manually builds taxonomy files without NCBI rsync dependency

set -euo pipefail

DB_PATH="/Users/kalinovdameus/Developer/Vibrion/data/kraken2_haiti_custom"
TAX_PATH="${DB_PATH}/taxonomy"

echo "Creating custom taxonomy for Kraken2..."
echo ""

# Create taxonomy directory
mkdir -p "${TAX_PATH}"

# Create nodes.dmp (taxonomy tree structure)
# Format: tax_id | parent_tax_id | rank | ... other fields
cat > "${TAX_PATH}/nodes.dmp" << 'EOF'
1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
2	|	131567	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
131567	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
1224	|	2	|	phylum	|	Pseudomonadota	|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
1236	|	1224	|	class	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
72274	|	1236	|	order	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
641	|	72274	|	family	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
662	|	641	|	genus	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
666	|	662	|	species	|	Vibrio cholerae	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
2010EL1786	|	666	|	strain	|	Vibrio cholerae 2010EL-1786 (Haiti)	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
N16961	|	666	|	strain	|	Vibrio cholerae N16961 (7th pandemic)	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
O139	|	666	|	strain	|	Vibrio cholerae O139	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
O37	|	666	|	strain	|	Vibrio cholerae O37	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
659	|	662	|	species	|	Vibrio mimicus	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
670	|	662	|	species	|	Vibrio parahaemolyticus	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
672	|	662	|	species	|	Vibrio vulnificus	|	0	|	1	|	11	|	0	|	0	|	0	|	0	|	0	|		|
EOF

# Create names.dmp (taxonomy names)
# Format: tax_id | name | unique name | name class
cat > "${TAX_PATH}/names.dmp" << 'EOF'
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
131567	|	cellular organisms	|		|	scientific name	|
1224	|	Pseudomonadota	|		|	scientific name	|
1224	|	Proteobacteria	|		|	synonym	|
1236	|	Gammaproteobacteria	|		|	scientific name	|
72274	|	Vibrionales	|		|	scientific name	|
641	|	Vibrionaceae	|		|	scientific name	|
662	|	Vibrio	|		|	scientific name	|
666	|	Vibrio cholerae	|		|	scientific name	|
2010EL1786	|	Vibrio cholerae 2010EL-1786	|		|	scientific name	|
2010EL1786	|	Haiti 2010 outbreak strain	|		|	genbank common name	|
N16961	|	Vibrio cholerae N16961	|		|	scientific name	|
N16961	|	7th pandemic reference strain	|		|	genbank common name	|
O139	|	Vibrio cholerae O139	|		|	scientific name	|
O37	|	Vibrio cholerae O37	|		|	scientific name	|
659	|	Vibrio mimicus	|		|	scientific name	|
670	|	Vibrio parahaemolyticus	|		|	scientific name	|
672	|	Vibrio vulnificus	|		|	scientific name	|
EOF

echo "✓ Created nodes.dmp and names.dmp"

# Create preliminary_accessions.txt for each library sequence
# This maps sequence IDs to taxonomy IDs
PRELIM_FILE="${DB_PATH}/library/added/prelim_map.txt"
mkdir -p "${DB_PATH}/library/added"

cat > "${PRELIM_FILE}" << 'EOF'
ACCESSION	TAXID
2010EL-1786	2010EL1786
CP003070.1	N16961
NC_022965.1	N16961
CP024868.1	O139
CP024869.1	O139
CP009262.1	O37
CP009263.1	O37
NZ_CP014039.1	659
NZ_CP014040.1	659
NC_004603.1	670
NC_004605.1	670
NC_004459.3	672
NC_005139.1	672
EOF

echo "✓ Created preliminary accessions map"

# Create seqid2taxid.map for the database
# This is what Kraken2 actually uses during build
SEQID_MAP="${TAX_PATH}/seqid2taxid.map"

cat > "${SEQID_MAP}" << 'EOF'
2010EL-1786	2010EL1786
CP003070.1	N16961
NC_022965.1	N16961
CP024868.1	O139
CP024869.1	O139
CP009262.1	O37
CP009263.1	O37
NZ_CP014039.1	659
NZ_CP014040.1	659
NC_004603.1	670
NC_004605.1	670
NC_004459.3	672
NC_005139.1	672
EOF

echo "✓ Created seqid2taxid.map"

# Update the library FASTA files to include NCBI-style sequence IDs
# Kraken2 expects headers like: >NC_012345.1 or >kraken:taxid|666|NC_012345.1
echo ""
echo "Updating FASTA headers with taxonomy..."

LIBRARY_PATH="/Users/kalinovdameus/Developer/Vibrion/data/kraken2_library"

# Process each genome and add to database library with proper headers
for category_dir in "${LIBRARY_PATH}"/*; do
    if [ -d "${category_dir}" ]; then
        category=$(basename "${category_dir}")
        echo "  Processing ${category}..."
        
        for fasta in "${category_dir}"/*.fasta; do
            if [ -f "${fasta}" ]; then
                genome=$(basename "${fasta}" .fasta)
                
                # Get taxid for this genome
                case "${genome}" in
                    "2010EL-1786") taxid="2010EL1786" ;;
                    "N16961") taxid="N16961" ;;
                    "O139") taxid="O139" ;;
                    "O37") taxid="O37" ;;
                    "mimicus") taxid="659" ;;
                    "parahaemolyticus") taxid="670" ;;
                    "vulnificus") taxid="672" ;;
                    *) taxid="666" ;;  # Generic V. cholerae
                esac
                
                # Create library entry with kraken-formatted headers
                out_fasta="${DB_PATH}/library/added/${genome}.fna"
                
                # Rewrite FASTA with kraken:taxid headers
                awk -v taxid="${taxid}" '
                    /^>/ { 
                        # Extract sequence ID (first word after >)
                        seqid = $1
                        gsub(/^>/, "", seqid)
                        # Create kraken-formatted header
                        print ">kraken:taxid|" taxid "|" seqid
                        next
                    }
                    { print }
                ' "${fasta}" > "${out_fasta}"
                
                echo "    Added ${genome} → taxid ${taxid}"
            fi
        done
    fi
done

echo ""
echo "✓ Taxonomy structure created successfully!"
echo ""
echo "Database ready for building with:"
echo "  kraken2-build --build --db ${DB_PATH}"
echo ""
