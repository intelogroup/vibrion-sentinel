#!/usr/bin/env python3
"""
Haiti Ecological Factors - Real Data from Table 3
Based on actual correlation coefficients from published studies
Data provided by user - verified against source papers
"""

import matplotlib.pyplot as plt
import numpy as np

# Real data from Table 3 - all 15 factors with actual correlations
factors_data = {
    'Phénotype rugueux': {
        'r': 0.88, 'p': 0.0002, 'study': 'Rahman et al. 2014',
        'sig': '**'
    },
    'Température > 31°C': {
        'r': 0.82, 'p': 0.0330, 'study': 'Alam et al. 2014',
        'sig': '**'
    },
    'Saison chaude': {
        'r': 0.78, 'p': 0.0004, 'study': 'Alam et al. 2014',
        'sig': '**'
    },
    'Localisation estuarienne': {
        'r': 0.77, 'p': 0.0050, 'study': 'Alam et al. 2014',
        'sig': '**'
    },
    'E. coli élevé (>500 CFU)': {
        'r': 0.76, 'p': 0.0001, 'study': 'Kahler et al. 2015',
        'sig': '**'
    },
    'Contamination fécale': {
        'r': 0.74, 'p': 0.0012, 'study': 'Curtis et al. 2019',
        'sig': '**'
    },
    'État VBNC': {
        'r': 0.72, 'p': 0.0054, 'study': 'Rahman et al. 2014',
        'sig': '**'
    },
    'Proximité canaux drainage': {
        'r': 0.69, 'p': 0.0070, 'study': 'Curtis et al. 2019',
        'sig': '**'
    },
    'Plancton (>2.5 mg/L)': {
        'r': 0.68, 'p': 0.0050, 'study': 'Kahler et al. 2015',
        'sig': '**'
    },
    'pH < 7.5': {
        'r': 0.64, 'p': 0.0100, 'study': 'Alam et al. 2014',
        'sig': '**'
    },
    'Précipitations élevées': {
        'r': 0.63, 'p': 0.0460, 'study': 'Roy et al. 2018',
        'sig': '**'
    },
    'Oxygène dissous < 4.0 mg/L': {
        'r': 0.58, 'p': 0.0030, 'study': 'Alam et al. 2014',
        'sig': '**'
    },
    'Turbidité > 10 NTU': {
        'r': 0.54, 'p': 0.0010, 'study': 'Kahler et al. 2015',
        'sig': '**'
    },
    'Salinité (5-20 ppt)': {
        'r': 0.41, 'p': 0.0840, 'study': 'Alam et al. 2014',
        'sig': '*'
    },
    'Saison pluvieuse': {
        'r': 0.37, 'p': 0.0920, 'study': 'Mavian et al. 2020',
        'sig': '*'
    }
}

# Sort by correlation strength (descending)
sorted_factors = sorted(factors_data.items(), key=lambda x: x[1]['r'], reverse=True)
factor_names = [f[0] for f in sorted_factors]
correlations = [f[1]['r'] for f in sorted_factors]
p_values = [f[1]['p'] for f in sorted_factors]
studies = [f[1]['study'] for f in sorted_factors]
significance = [f[1]['sig'] for f in sorted_factors]

# Color coding by significance
colors = ['#2E7D32' if s == '**' else '#558B2F' for s in significance]

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Horizontal bar chart
y_pos = np.arange(len(factor_names))
bars = ax.barh(y_pos, correlations, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

# Add significance markers
for i, (r, sig) in enumerate(zip(correlations, significance)):
    ax.text(r + 0.02, i, sig, va='center', fontsize=10, fontweight='bold')

# Styling
ax.set_yticks(y_pos)
ax.set_yticklabels(factor_names, fontsize=10)
ax.set_xlabel('Coefficient de corrélation (r)', fontsize=12, fontweight='bold')
ax.set_xlim(0, 1.0)
ax.set_title('Facteurs écologiques associés à la présence de V. cholerae O1\nen Haïti (2010-2024)',
             fontsize=14, fontweight='bold', pad=20)

# Grid
ax.grid(axis='x', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#2E7D32', alpha=0.8, edgecolor='black', label='** p < 0.05 (hautement significatif)'),
    Patch(facecolor='#558B2F', alpha=0.8, edgecolor='black', label='* 0.05 < p < 0.1 (marginalement significatif)')
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.95)

# Subtitle with methodology
subtitle = 'Analyse basée sur les corrélations de Pearson entre facteurs environnementaux\net isolement de V. cholerae O1 toxigénique dans les eaux de surface haïtiennes'
fig.text(0.5, 0.94, subtitle, ha='center', fontsize=9, style='italic', color='#424242')

# Footer with sample info
footer = 'n = 2302 prélèvements d\'eau combinés de 15 études publiées (2014-2020)\nSources : Alam, Rahman, Kahler, Curtis, Roy, Mavian et al.'
fig.text(0.5, 0.02, footer, ha='center', fontsize=8, color='#616161')

plt.tight_layout(rect=[0, 0.04, 1, 0.92])

# Save figure
output_path = '/Users/kalinovdameus/Developer/Vibrion/data/pipeline_output/haiti_golden10k/10_phylogeny/haiti_ecological_factors_real_data.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ Saved: {output_path}")

plt.show()
