"""
REA — Restriction Enzyme Analyzer
A professional Biopython-powered tool for restriction mapping,
cut site detection, and virtual gel electrophoresis simulation.
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import (
    RestrictionBatch, AllEnzymes, CommOnly,
    EcoRI, HindIII, BamHI, NotI, XhoI,
    NcoI, SalI, XbaI, PstI, SphI, KpnI, SacI,
    EcoRV, ClaI, NheI, SpeI, ApaI, MluI
)
import io
import math
import re

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="REA · Restriction Enzyme Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@300;400;500&display=swap');

:root {
    --bg-primary: #0B1220;
    --bg-secondary: #111827;
    --bg-card: #151d2e;
    --bg-elevated: #1a2332;
    --accent-blue: #4F7DF3;
    --accent-green: #22C55E;
    --accent-teal: #14B8A6;
    --text-primary: #F8FAFC;
    --text-secondary: #CBD5E1;
    --text-muted: #94A3B8;
    --border-color: #1e293b;
    --border-subtle: rgba(255,255,255,0.06);
}

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    font-weight: 400;
}

/* Global dark theme */
.stApp {
    background: linear-gradient(180deg, #0B1220 0%, #0f172a 100%);
}

/* Hide default streamlit chrome but keep sidebar toggle */
#MainMenu, footer { visibility: hidden; }
header[data-testid="stHeader"] { background: transparent; }
button[data-testid="stBaseButton-headerNoPadding"] { visibility: visible !important; }
.block-container { 
    padding-top: 1rem; 
    padding-bottom: 2rem;
    max-width: 1400px;
}

/* Topbar */
.topbar {
    background: linear-gradient(135deg, rgba(17,24,39,0.95) 0%, rgba(11,18,32,0.98) 100%);
    backdrop-filter: blur(10px);
    color: var(--text-primary);
    padding: 12px 24px;
    margin: -1rem -1rem 1.5rem -1rem;
    display: flex;
    align-items: center;
    justify-content: space-between;
    border-bottom: 1px solid var(--border-color);
    box-shadow: 0 4px 20px rgba(0,0,0,0.3);
}
.topbar-brand {
    font-family: 'Inter', sans-serif;
    font-size: 14px;
    font-weight: 600;
    letter-spacing: 0.02em;
    color: var(--text-primary);
}
.topbar-brand span {
    color: var(--accent-blue);
}
.topbar-sub {
    font-family: 'JetBrains Mono', monospace;
    font-size: 10px;
    color: var(--text-muted);
    letter-spacing: 0.03em;
    margin-top: 2px;
}

/* Section headers */
.section-label {
    font-family: 'Inter', sans-serif;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: var(--text-muted);
    margin-bottom: 12px;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--border-subtle);
    display: flex;
    align-items: center;
    gap: 8px;
}
.section-label::before {
    content: '';
    width: 3px;
    height: 12px;
    background: var(--accent-blue);
    border-radius: 2px;
}

/* Metric cards */
.metric-row {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 16px;
    margin: 20px 0;
}
.metric-card {
    background: linear-gradient(145deg, rgba(21,29,46,0.8) 0%, rgba(17,24,39,0.9) 100%);
    border: 1px solid var(--border-color);
    padding: 20px 22px;
    border-radius: 12px;
    position: relative;
    overflow: hidden;
    box-shadow: 0 4px 16px rgba(0,0,0,0.2), inset 0 1px 0 rgba(255,255,255,0.03);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
}
.metric-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 24px rgba(0,0,0,0.3), inset 0 1px 0 rgba(255,255,255,0.05);
}
.metric-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 3px;
    background: linear-gradient(90deg, var(--accent-blue), var(--accent-teal));
    opacity: 0.8;
}
.metric-card .val {
    font-family: 'Inter', sans-serif;
    font-size: 28px;
    font-weight: 700;
    color: var(--text-primary);
    line-height: 1.1;
    letter-spacing: -0.02em;
}
.metric-card .lbl {
    font-size: 12px;
    color: var(--text-muted);
    margin-top: 6px;
    letter-spacing: 0.01em;
    font-weight: 500;
}

/* Enzyme badges */
.enzyme-badge {
    display: inline-flex;
    align-items: center;
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    font-weight: 500;
    padding: 6px 12px;
    border-radius: 6px;
    margin: 3px;
    transition: all 0.15s ease;
    cursor: default;
}
.enzyme-badge.cutter {
    background: linear-gradient(135deg, rgba(34,197,94,0.15) 0%, rgba(20,184,166,0.1) 100%);
    color: #4ade80;
    border: 1px solid rgba(34,197,94,0.3);
    box-shadow: 0 2px 8px rgba(34,197,94,0.1);
}
.enzyme-badge.cutter:hover {
    background: linear-gradient(135deg, rgba(34,197,94,0.25) 0%, rgba(20,184,166,0.15) 100%);
    transform: translateY(-1px);
}
.enzyme-badge.cutter::before {
    content: '✓';
    margin-right: 6px;
    font-size: 10px;
}
.enzyme-badge.nocutter {
    background: rgba(148,163,184,0.08);
    color: var(--text-muted);
    border: 1px solid rgba(148,163,184,0.15);
}
.enzyme-badge.nocutter:hover {
    background: rgba(148,163,184,0.12);
}

/* Sequence display */
.seq-block {
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
    font-weight: 400;
    background: linear-gradient(145deg, #0a0f1a 0%, #0d1320 100%);
    color: #6ee7b7;
    padding: 20px;
    border-radius: 10px;
    border: 1px solid var(--border-color);
    overflow-x: auto;
    line-height: 1.9;
    letter-spacing: 0.08em;
    white-space: pre-wrap;
    word-break: break-all;
    max-height: 220px;
    overflow-y: auto;
    box-shadow: inset 0 2px 10px rgba(0,0,0,0.3);
}

/* Info box */
.info-box {
    background: linear-gradient(135deg, rgba(34,197,94,0.08) 0%, rgba(20,184,166,0.05) 100%);
    border-left: 3px solid var(--accent-green);
    border-radius: 0 8px 8px 0;
    padding: 12px 16px;
    font-size: 13px;
    color: #86efac;
    margin: 10px 0;
    border: 1px solid rgba(34,197,94,0.15);
    border-left: 3px solid var(--accent-green);
}
.warn-box {
    background: linear-gradient(135deg, rgba(245,158,11,0.08) 0%, rgba(239,68,68,0.05) 100%);
    border-left: 3px solid #f59e0b;
    border-radius: 0 8px 8px 0;
    padding: 12px 16px;
    font-size: 13px;
    color: #fcd34d;
    margin: 10px 0;
    border: 1px solid rgba(245,158,11,0.15);
    border-left: 3px solid #f59e0b;
}

/* Sidebar */
section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0d1525 0%, #0a1018 100%);
    border-right: 1px solid var(--border-color);
}
section[data-testid="stSidebar"] > div {
    padding-top: 1rem;
}
section[data-testid="stSidebar"] * { 
    color: var(--text-secondary) !important; 
}
section[data-testid="stSidebar"] .stSelectbox label,
section[data-testid="stSidebar"] .stTextInput label,
section[data-testid="stSidebar"] .stTextArea label,
section[data-testid="stSidebar"] .stFileUploader label,
section[data-testid="stSidebar"] .stMultiSelect label,
section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 {
    font-family: 'Inter', sans-serif !important;
    font-size: 11px !important;
    font-weight: 600 !important;
    letter-spacing: 0.06em !important;
    text-transform: uppercase !important;
    color: var(--text-muted) !important;
}

/* Sidebar inputs */
section[data-testid="stSidebar"] .stSelectbox > div > div,
section[data-testid="stSidebar"] .stTextInput > div > div > input,
section[data-testid="stSidebar"] .stTextArea > div > div > textarea {
    background: var(--bg-elevated) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 8px !important;
    color: var(--text-primary) !important;
}
section[data-testid="stSidebar"] .stSelectbox > div > div:hover,
section[data-testid="stSidebar"] .stTextInput > div > div > input:hover {
    border-color: var(--accent-blue) !important;
}

/* Sidebar button */
section[data-testid="stSidebar"] .stButton > button {
    background: linear-gradient(135deg, var(--accent-blue) 0%, #3b6ce7 100%) !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    font-family: 'Inter', sans-serif !important;
    font-size: 12px !important;
    font-weight: 600 !important;
    letter-spacing: 0.03em !important;
    text-transform: uppercase !important;
    width: 100% !important;
    padding: 12px 16px !important;
    margin-top: 12px !important;
    box-shadow: 0 4px 14px rgba(79,125,243,0.3) !important;
    transition: all 0.2s ease !important;
}
section[data-testid="stSidebar"] .stButton > button:hover {
    background: linear-gradient(135deg, #5b8af5 0%, #4a7cf0 100%) !important;
    transform: translateY(-1px) !important;
    box-shadow: 0 6px 20px rgba(79,125,243,0.4) !important;
}
section[data-testid="stSidebar"] .stButton > button:active {
    transform: scale(0.98) !important;
}

/* Multi-select tags */
section[data-testid="stSidebar"] .stMultiSelect span[data-baseweb="tag"] {
    background: linear-gradient(135deg, rgba(79,125,243,0.2) 0%, rgba(59,108,231,0.15) 100%) !important;
    border: 1px solid rgba(79,125,243,0.3) !important;
    border-radius: 6px !important;
    color: #93c5fd !important;
}

/* Tabs */
.stTabs [data-baseweb="tab-list"] {
    gap: 0;
    background: transparent;
    border-bottom: 1px solid var(--border-color);
    padding: 0 4px;
}
.stTabs [data-baseweb="tab"] {
    font-family: 'Inter', sans-serif;
    font-size: 12px;
    font-weight: 500;
    letter-spacing: 0.02em;
    padding: 12px 20px;
    border-radius: 8px 8px 0 0;
    border-bottom: 2px solid transparent;
    color: var(--text-muted);
    background: transparent;
    transition: all 0.15s ease;
}
.stTabs [data-baseweb="tab"]:hover {
    color: var(--text-secondary);
    background: rgba(79,125,243,0.05);
}
.stTabs [aria-selected="true"] {
    color: var(--text-primary) !important;
    border-bottom-color: var(--accent-blue) !important;
    background: rgba(79,125,243,0.08) !important;
}

/* Dataframe */
.stDataFrame { 
    border: 1px solid var(--border-color) !important; 
    border-radius: 10px !important;
    overflow: hidden;
}
.stDataFrame [data-testid="stDataFrameResizable"] {
    background: var(--bg-card) !important;
}

/* Slider */
section[data-testid="stSidebar"] .stSlider > div > div > div {
    background: var(--border-color) !important;
}
section[data-testid="stSidebar"] .stSlider > div > div > div > div {
    background: var(--accent-blue) !important;
}

/* Dividers */
hr { 
    border: none; 
    border-top: 1px solid var(--border-color); 
    margin: 24px 0; 
    opacity: 0.5;
}

/* Landing page */
.landing-container {
    text-align: center;
    padding: 60px 20px;
}
.landing-title {
    font-family: 'Inter', sans-serif;
    font-size: 13px;
    font-weight: 600;
    letter-spacing: 0.15em;
    text-transform: uppercase;
    color: var(--accent-blue);
    margin-bottom: 20px;
}
.landing-desc {
    font-size: 18px;
    color: var(--text-primary);
    font-weight: 400;
    line-height: 1.7;
    margin-bottom: 28px;
    max-width: 500px;
    margin-left: auto;
    margin-right: auto;
}
.landing-features {
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
    color: var(--text-muted);
    line-height: 2.2;
}
.landing-features span {
    color: var(--accent-teal);
}

/* Plotly chart container */
.js-plotly-plot {
    border-radius: 12px;
    overflow: hidden;
}

/* Scrollbar styling */
::-webkit-scrollbar {
    width: 8px;
    height: 8px;
}
::-webkit-scrollbar-track {
    background: var(--bg-primary);
}
::-webkit-scrollbar-thumb {
    background: var(--border-color);
    border-radius: 4px;
}
::-webkit-scrollbar-thumb:hover {
    background: #334155;
}

/* File uploader */
section[data-testid="stSidebar"] .stFileUploader > div {
    background: var(--bg-elevated) !important;
    border: 1px dashed var(--border-color) !important;
    border-radius: 8px !important;
}
section[data-testid="stSidebar"] .stFileUploader > div:hover {
    border-color: var(--accent-blue) !important;
}
</style>
""", unsafe_allow_html=True)

# ── Constants ─────────────────────────────────────────────────────────────────
COMMON_ENZYMES = {
    "EcoRI":  EcoRI,  "HindIII": HindIII, "BamHI":  BamHI,
    "NotI":   NotI,   "XhoI":    XhoI,    "NcoI":   NcoI,
    "SalI":   SalI,   "XbaI":    XbaI,    "PstI":   PstI,
    "SphI":   SphI,   "KpnI":    KpnI,    "SacI":   SacI,
    "EcoRV":  EcoRV,  "ClaI":    ClaI,    "NheI":   NheI,
    "SpeI":   SpeI,   "ApaI":    ApaI,    "MluI":   MluI,
}

SAMPLE_SEQUENCES = {
    "pUC19 MCS region (synthetic)": (
        "AATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAG"
    ),
    "Lambda phage fragment (synthetic)": (
        "GGATCCAAGCTTAAGCTTGAATTCCTCGAGGTCGACGCATGCGAGCTCAAGCTTGAATTCGGATCCCTCGAGGTCGACGCATGCGAATTCAAGCTTGGATCCCTCGAGAAGCTTGAATTCGGATCCGTCGACGCATGCGAGCTCAAGCTT"
    ),
    "Custom BRCA1 exon region (synthetic)": (
        "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAG"
        "AGAGTCCCTGTGGATTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGAGTCC"
        "CTGTGGAATTCAAGCTTGGATCCGAATTCAAGCTTGGATCC"
    ),
}

DNA_LADDER = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
              1200, 1500, 2000, 2500, 3000]

# ── Helpers ───────────────────────────────────────────────────────────────────
def validate_dna(seq: str) -> tuple[bool, str]:
    seq = seq.upper().strip()
    invalid = set(seq) - set("ATGCNRYWSMKHBVD\n ")
    if invalid:
        return False, f"Invalid characters: {', '.join(sorted(invalid))}"
    clean = re.sub(r'[\s\n]', '', seq)
    if len(clean) < 20:
        return False, "Sequence must be at least 20 bp."
    return True, clean


def analyze_sequence(seq_str: str, selected_enzyme_names: list) -> dict:
    seq = Seq(seq_str)
    enzymes = [COMMON_ENZYMES[e] for e in selected_enzyme_names]
    rb = RestrictionBatch(enzymes)
    results = rb.search(seq, linear=True)

    analysis = {}
    for enz_obj, positions in results.items():
        name = str(enz_obj)
        cuts = sorted(positions)
        fragments = _compute_fragments(len(seq_str), cuts)
        analysis[name] = {
            "cuts": cuts,
            "cut_count": len(cuts),
            "fragments": fragments,
            "recognition_site": enz_obj.site,
            "overhang": _get_overhang(enz_obj),
            "cuts_sequence": len(cuts) > 0,
        }
    return analysis


def _compute_fragments(length: int, cuts: list) -> list:
    if not cuts:
        return [length]
    boundaries = [0] + cuts + [length]
    return sorted([boundaries[i+1] - boundaries[i]
                   for i in range(len(boundaries)-1)], reverse=True)


def _get_overhang(enz) -> str:
    try:
        overhang = enz.ovhg
        if overhang > 0:
            return f"5' overhang ({overhang} nt)"
        elif overhang < 0:
            return f"3' overhang ({abs(overhang)} nt)"
        else:
            return "Blunt end"
    except:
        return "Unknown"


def gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq) * 100, 2) if seq else 0.0


def fetch_from_ncbi(accession: str, email: str) -> SeqRecord | None:
    Entrez.email = email
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession,
            rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        st.error(f"NCBI fetch failed: {e}")
        return None


def format_sequence_display(seq: str, width: int = 60) -> str:
    lines = []
    for i in range(0, len(seq), width):
        chunk = seq[i:i+width]
        pos = str(i+1).rjust(6)
        lines.append(f"{pos}  {chunk}")
    return "\n".join(lines)


# ── Plotly Charts ─────────────────────────────────────────────────────────────

PLOT_LAYOUT = dict(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='#111827',
    font=dict(family='Inter, sans-serif', size=11, color='#CBD5E1'),
    margin=dict(l=50, r=30, t=50, b=50),
    xaxis=dict(gridcolor='#1e293b', linecolor='#334155', tickfont=dict(size=10, color='#94A3B8')),
    yaxis=dict(gridcolor='#1e293b', linecolor='#334155', tickfont=dict(size=10, color='#94A3B8')),
)


def plot_cut_map(seq_len: int, analysis: dict, seq_str: str):
    enzymes_with_cuts = {k: v for k, v in analysis.items() if v['cut_count'] > 0}
    if not enzymes_with_cuts:
        st.info("No enzymes cut this sequence.")
        return

    fig = go.Figure()
    colors = ['#4F7DF3', '#22C55E', '#14B8A6', '#F59E0B', '#EC4899',
              '#8B5CF6', '#06B6D4', '#84CC16', '#F97316', '#6366F1']

    # Backbone
    fig.add_shape(type="line", x0=0, x1=seq_len, y0=0, y1=0,
                  line=dict(color='#64748b', width=3))

    # Add start/end markers
    fig.add_annotation(x=0, y=0, text="5'", showarrow=False,
                       font=dict(size=10, color='#94A3B8'), xanchor='right', xshift=-6)
    fig.add_annotation(x=seq_len, y=0, text="3'", showarrow=False,
                       font=dict(size=10, color='#94A3B8'), xanchor='left', xshift=6)

    y_levels = list(range(1, len(enzymes_with_cuts) + 1))
    for idx, (name, data) in enumerate(enzymes_with_cuts.items()):
        color = colors[idx % len(colors)]
        y = y_levels[idx]
        for pos in data['cuts']:
            # Vertical cut tick
            fig.add_shape(type="line", x0=pos, x1=pos, y0=-0.3, y1=0.3,
                          line=dict(color=color, width=2))
            # Leader to label row
            fig.add_shape(type="line", x0=pos, x1=pos, y0=0.3, y1=y - 0.15,
                          line=dict(color=color, width=0.8, dash='dot'))
            fig.add_trace(go.Scatter(
                x=[pos], y=[y], mode='markers+text',
                marker=dict(size=8, color=color, symbol='diamond'),
                text=[f"{pos}"],
                textposition='top center',
                textfont=dict(size=9, color=color),
                name=name, showlegend=False,
                hovertemplate=f"<b>{name}</b><br>Position: {pos} bp<extra></extra>"
            ))
        # Enzyme label on left
        fig.add_annotation(
            x=0, y=y, text=name, showarrow=False,
            font=dict(size=10, color=color, family='Inter'),
            xanchor='right', xshift=-10
        )

    max_y = len(enzymes_with_cuts) + 0.8
    fig.update_layout(
        **{k: v for k, v in PLOT_LAYOUT.items() if k not in ('xaxis', 'yaxis')},
        height=max(280, len(enzymes_with_cuts) * 55 + 100),
        title=dict(text=f"Cut Site Map · {seq_len} bp", font=dict(size=13, color='#F8FAFC'), x=0),
        xaxis=dict(
            title="Position (bp)",
            range=[-seq_len * 0.12, seq_len * 1.05],
            **PLOT_LAYOUT['xaxis']
        ),
        yaxis=dict(
            visible=False,
            range=[-0.8, max_y],
        ),
        showlegend=False,
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_gel(analysis: dict, seq_len: int):
    enzymes_with_cuts = {k: v for k, v in analysis.items() if v['cut_count'] > 0}
    if not enzymes_with_cuts:
        st.info("No fragments to display — no enzymes cut this sequence.")
        return

    all_lane_names = ['Ladder'] + list(enzymes_with_cuts.keys())
    num_lanes = len(all_lane_names)

    fig = go.Figure()

    lane_width = 0.6
    gel_color = '#0a0f1a'
    band_colors = ['#22C55E', '#14B8A6', '#4F7DF3', '#F59E0B', 
                   '#EC4899', '#8B5CF6', '#06B6D4', '#84CC16']
    # Glow colors (rgba with low opacity)
    glow_colors = ['rgba(34,197,94,0.15)', 'rgba(20,184,166,0.15)', 'rgba(79,125,243,0.15)', 
                   'rgba(245,158,11,0.15)', 'rgba(236,72,153,0.15)', 'rgba(139,92,246,0.15)',
                   'rgba(6,182,212,0.15)', 'rgba(132,204,22,0.15)']

    def bp_to_y(bp):
        if bp <= 0:
            return 0
        return math.log10(bp)

    # Gel background with subtle gradient effect
    y_min = bp_to_y(50)
    y_max = bp_to_y(4000)
    fig.add_shape(type="rect", x0=0.1, x1=num_lanes + 0.9,
                  y0=y_min - 0.05, y1=y_max + 0.15,
                  fillcolor=gel_color, 
                  line=dict(color='#1e293b', width=1))

    # Lane separators
    for i in range(1, num_lanes):
        fig.add_shape(type="line", x0=i + 0.65, x1=i + 0.65,
                      y0=y_min - 0.05, y1=y_max + 0.15,
                      line=dict(color='#1e293b', width=0.5))

    # Ladder lane with glow effect
    for bp in DNA_LADDER:
        y = bp_to_y(bp)
        if y < y_min or y > y_max + 0.1:
            continue
        intensity = 0.7 + 0.2 * (1 - (bp / max(DNA_LADDER)))
        # Soft glow behind band
        fig.add_shape(type="rect",
                      x0=1 - lane_width / 2 - 0.02, x1=1 + lane_width / 2 + 0.02,
                      y0=y - 0.018, y1=y + 0.018,
                      fillcolor=f'rgba(148,163,184,0.15)',
                      line_width=0)
        # Main band
        fig.add_shape(type="rect",
                      x0=1 - lane_width / 2, x1=1 + lane_width / 2,
                      y0=y - 0.012, y1=y + 0.012,
                      fillcolor=f'rgba(226,232,240,{intensity:.2f})',
                      line_width=0)
        fig.add_annotation(x=0.38, y=y, text=f"{bp}",
                           showarrow=False,
                           font=dict(size=8, color='#94A3B8',
                                     family='Inter'),
                           xanchor='right')

    # Sample lanes with glow
    for lane_idx, (enz_name, data) in enumerate(enzymes_with_cuts.items()):
        lane_x = lane_idx + 2
        color = band_colors[lane_idx % len(band_colors)]
        glow_color = glow_colors[lane_idx % len(glow_colors)]
        fragments = data['fragments']
        for frag_bp in fragments:
            if frag_bp < 50:
                continue
            y = bp_to_y(frag_bp)
            if y < y_min - 0.02 or y > y_max + 0.12:
                continue
            thickness = max(0.008, 0.028 - frag_bp / 120000)
            # Soft glow behind band
            fig.add_shape(type="rect",
                          x0=lane_x - lane_width / 2 - 0.02,
                          x1=lane_x + lane_width / 2 + 0.02,
                          y0=y - thickness - 0.006, y1=y + thickness + 0.006,
                          fillcolor=glow_color,
                          line_width=0)
            # Main band
            fig.add_shape(type="rect",
                          x0=lane_x - lane_width / 2,
                          x1=lane_x + lane_width / 2,
                          y0=y - thickness, y1=y + thickness,
                          fillcolor=color,
                          line_width=0)
            fig.add_trace(go.Scatter(
                x=[lane_x], y=[y], mode='markers',
                marker=dict(size=0.1, opacity=0),
                hovertemplate=f"<b>{enz_name}</b><br>{frag_bp} bp<extra></extra>",
                showlegend=False
            ))

    # Lane labels
    tick_vals = list(range(1, num_lanes + 1))
    tick_texts = all_lane_names

    ytick_vals = [bp_to_y(bp) for bp in [100, 200, 500, 1000, 2000, 3000]
                  if y_min <= bp_to_y(bp) <= y_max + 0.1]
    ytick_texts = ['100', '200', '500', '1000', '2000', '3000']
    ytick_texts = [t for t, v in zip(ytick_texts, [bp_to_y(bp) for bp in [100, 200, 500, 1000, 2000, 3000]])
                   if y_min <= v <= y_max + 0.1]

    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(family='Inter, sans-serif', size=10, color='#CBD5E1'),
        height=500,
        margin=dict(l=65, r=30, t=55, b=65),
        title=dict(text="Virtual Gel Electrophoresis", font=dict(size=13, color='#F8FAFC'), x=0),
        xaxis=dict(
            tickvals=tick_vals, ticktext=tick_texts,
            tickfont=dict(size=10, color='#94A3B8'),
            gridcolor='rgba(0,0,0,0)', linecolor='rgba(0,0,0,0)',
            range=[0.1, num_lanes + 0.9],
        ),
        yaxis=dict(
            title=dict(text="Fragment Size (bp)", font=dict(size=11, color='#94A3B8')),
            tickvals=ytick_vals, ticktext=ytick_texts,
            tickfont=dict(size=10, color='#94A3B8'),
            gridcolor='rgba(0,0,0,0)', linecolor='rgba(0,0,0,0)',
            range=[y_min - 0.05, y_max + 0.15],
        ),
        showlegend=False,
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_gc_window(seq: str, window: int = 50):
    gc_vals, positions = [], []
    for i in range(0, len(seq) - window, max(1, window // 4)):
        chunk = seq[i:i+window]
        gc_vals.append(gc_content(chunk))
        positions.append(i + window // 2)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=positions, y=gc_vals, mode='lines',
        line=dict(color='#22C55E', width=2, shape='spline'),
        fill='tozeroy', fillcolor='rgba(34,197,94,0.1)',
        hovertemplate="Position: %{x} bp<br>GC: %{y:.1f}%<extra></extra>"
    ))
    fig.add_hline(y=50, line=dict(color='#64748b', width=1, dash='dash'),
                  annotation_text="50%", annotation_font_size=10, annotation_font_color='#94A3B8')
    fig.update_layout(
        **{k: v for k, v in PLOT_LAYOUT.items() if k not in ('xaxis', 'yaxis')},
        height=240,
        title=dict(text=f"GC Content · {window} bp sliding window", font=dict(size=13, color='#F8FAFC'), x=0),
        xaxis=dict(title="Position (bp)", **PLOT_LAYOUT['xaxis']),
        yaxis=dict(title="GC %", range=[0, 100], **PLOT_LAYOUT['yaxis']),
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_fragment_distribution(analysis: dict):
    cutters = {k: v for k, v in analysis.items() if v['cut_count'] > 0}
    if not cutters:
        return

    fig = go.Figure()
    colors = ['#4F7DF3', '#22C55E', '#14B8A6', '#F59E0B',
              '#EC4899', '#8B5CF6', '#06B6D4', '#84CC16']

    for idx, (name, data) in enumerate(cutters.items()):
        frags = data['fragments']
        color = colors[idx % len(colors)]
        fig.add_trace(go.Bar(
            name=name,
            x=[f"F{i+1}" for i in range(len(frags))],
            y=frags,
            marker=dict(color=color, line=dict(width=0)),
            hovertemplate=f"<b>{name}</b><br>Fragment: %{{x}}<br>Size: %{{y}} bp<extra></extra>"
        ))

    fig.update_layout(
        **{k: v for k, v in PLOT_LAYOUT.items() if k not in ('xaxis', 'yaxis')},
        height=300,
        barmode='group',
        title=dict(text="Fragment Size Distribution", font=dict(size=13, color='#F8FAFC'), x=0),
        xaxis=dict(title="Fragments (sorted by size)", **PLOT_LAYOUT['xaxis']),
        yaxis=dict(title="Size (bp)", **PLOT_LAYOUT['yaxis']),
        legend=dict(font=dict(size=10, color='#CBD5E1'), bgcolor='rgba(17,24,39,0.8)',
                    bordercolor='#1e293b', borderwidth=1)
    )
    st.plotly_chart(fig, use_container_width=True)


# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown('<p class="section-label">Input Source</p>', unsafe_allow_html=True)
    input_mode = st.selectbox(
        "input_mode", ["Manual / Paste", "Upload File", "Sample Sequences", "NCBI Accession"],
        label_visibility="collapsed"
    )
    st.markdown("---")

    seq_str = None

    if input_mode == "Manual / Paste":
        st.markdown('<p class="section-label">DNA Sequence</p>', unsafe_allow_html=True)
        raw_seq = st.text_area(
            "seq", height=160, label_visibility="collapsed",
            placeholder="Paste raw DNA sequence (ATGC...).\nSpaces and newlines are stripped automatically."
        )
        if raw_seq.strip():
            ok, result = validate_dna(raw_seq)
            if ok:
                seq_str = result
            else:
                st.markdown(f'<div class="warn-box">⚠ {result}</div>', unsafe_allow_html=True)

    elif input_mode == "Upload File":
        st.markdown('<p class="section-label">Upload FASTA/GenBank/Text</p>', unsafe_allow_html=True)
        uploaded_file = st.file_uploader(
            "file", type=["fasta", "fa", "fna", "gb", "gbk", "genbank", "txt"],
            label_visibility="collapsed"
        )
        if uploaded_file is not None:
            file_content = uploaded_file.read().decode("utf-8")
            file_name = uploaded_file.name.lower()

            try:
                from io import StringIO
                if file_name.endswith((".gb", ".gbk", ".genbank")):
                    records = list(SeqIO.parse(StringIO(file_content), "genbank"))
                    if not records:
                        raise ValueError("No sequences found in file")
                    record = records[0]
                    raw_seq = str(record.seq)
                    desc = record.description[:40] + "…" if len(record.description) > 40 else record.description
                    info = f"GenBank: {desc}"
                    if len(records) > 1:
                        info += f" (using 1 of {len(records)} records)"
                    st.markdown(f'<div class="info-box">{info}</div>', unsafe_allow_html=True)
                elif file_name.endswith((".fasta", ".fa", ".fna")):
                    records = list(SeqIO.parse(StringIO(file_content), "fasta"))
                    if not records:
                        raise ValueError("No sequences found in file")
                    record = records[0]
                    raw_seq = str(record.seq)
                    desc = record.description[:40] + "…" if len(record.description) > 40 else record.description
                    info = f"FASTA: {desc}"
                    if len(records) > 1:
                        info += f" (using 1 of {len(records)} records)"
                    st.markdown(f'<div class="info-box">{info}</div>', unsafe_allow_html=True)
                else:
                    raw_seq = file_content

                ok, result = validate_dna(raw_seq)
                if ok:
                    seq_str = result
                    st.markdown(
                        f'<div class="info-box">Loaded · {len(seq_str)} bp</div>',
                        unsafe_allow_html=True
                    )
                else:
                    st.markdown(f'<div class="warn-box">⚠ {result}</div>', unsafe_allow_html=True)
            except Exception as e:
                st.markdown(f'<div class="warn-box">⚠ Could not parse file: {e}</div>', unsafe_allow_html=True)

    elif input_mode == "Sample Sequences":
        st.markdown('<p class="section-label">Sample</p>', unsafe_allow_html=True)
        sample_name = st.selectbox("sample", list(SAMPLE_SEQUENCES.keys()),
                                   label_visibility="collapsed")
        seq_str = SAMPLE_SEQUENCES[sample_name]
        st.markdown(
            f'<div class="info-box">Sample loaded · {len(seq_str)} bp</div>',
            unsafe_allow_html=True
        )

    else:  # NCBI
        st.markdown('<p class="section-label">Accession Number</p>', unsafe_allow_html=True)
        accession = st.text_input("acc", placeholder="e.g. NM_007294",
                                  label_visibility="collapsed")
        st.markdown('<p class="section-label">Your Email (NCBI required)</p>',
                    unsafe_allow_html=True)
        email = st.text_input("email", placeholder="your@email.com",
                              label_visibility="collapsed")
        if st.button("Fetch from NCBI"):
            if accession and email:
                with st.spinner("Fetching..."):
                    record = fetch_from_ncbi(accession, email)
                    if record:
                        seq_str = str(record.seq)
                        st.session_state['ncbi_record'] = record
                        st.session_state['ncbi_seq'] = seq_str
                        st.markdown(
                            f'<div class="info-box">Fetched: {record.description[:50]}… · {len(seq_str)} bp</div>',
                            unsafe_allow_html=True
                        )
            else:
                st.warning("Enter accession and email.")

        if 'ncbi_seq' in st.session_state:
            seq_str = st.session_state['ncbi_seq']

    st.markdown("---")
    st.markdown('<p class="section-label">Enzyme Panel</p>', unsafe_allow_html=True)
    selected_enzymes = st.multiselect(
        "enzymes", list(COMMON_ENZYMES.keys()),
        default=["EcoRI", "HindIII", "BamHI", "NotI", "XhoI"],
        label_visibility="collapsed"
    )

    st.markdown("---")
    st.markdown('<p class="section-label">GC Window (bp)</p>', unsafe_allow_html=True)
    gc_window = st.slider("gcw", 20, 200, 50, 10, label_visibility="collapsed")

    st.markdown("---")
    run = st.button("Run Analysis", use_container_width=True)

# ── Main panel ────────────────────────────────────────────────────────────────
st.markdown("""
<div class="topbar">
  <div>
    <div class="topbar-brand"><span>REA</span> · Restriction Enzyme Analyzer</div>
    <div class="topbar-sub">Powered by Biopython · Entrez · Plotly</div>
  </div>
</div>
""", unsafe_allow_html=True)

if not seq_str:
    # Landing state
    c1, c2, c3 = st.columns([1, 2, 1])
    with c2:
        st.markdown("---")
        st.markdown("""
        <div class="landing-container">
          <div class="landing-title">Restriction Enzyme Analyzer</div>
          <div class="landing-desc">
            Paste a DNA sequence or fetch from NCBI,<br>
            select your enzyme panel, and run analysis.
          </div>
          <div class="landing-features">
            <span>◆</span> Cut Site Mapping &nbsp;&nbsp; <span>◆</span> Virtual Gel Electrophoresis<br>
            <span>◆</span> Fragment Analysis &nbsp;&nbsp; <span>◆</span> GC Content Profiling<br>
            <span>◆</span> Enzyme Overhang Classification
          </div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("---")

elif not run and seq_str:
    gcv = gc_content(seq_str)
    st.markdown(f"""
    <div class="metric-row">
      <div class="metric-card">
        <div class="val">{len(seq_str):,}</div>
        <div class="lbl">Sequence length (bp)</div>
      </div>
      <div class="metric-card">
        <div class="val">{gcv}%</div>
        <div class="lbl">GC content</div>
      </div>
      <div class="metric-card">
        <div class="val">{len(selected_enzymes)}</div>
        <div class="lbl">Enzymes selected</div>
      </div>
      <div class="metric-card">
        <div class="val">—</div>
        <div class="lbl">Click Run Analysis</div>
      </div>
    </div>
    """, unsafe_allow_html=True)
    st.markdown('<p class="section-label">Sequence Preview</p>', unsafe_allow_html=True)
    preview = format_sequence_display(seq_str[:300])
    st.markdown(f'<div class="seq-block">{preview}{"…" if len(seq_str) > 300 else ""}</div>',
                unsafe_allow_html=True)
    plot_gc_window(seq_str, gc_window)

else:
    if not selected_enzymes:
        st.warning("Select at least one enzyme from the sidebar.")
        st.stop()

    with st.spinner("Analysing sequence…"):
        analysis = analyze_sequence(seq_str, selected_enzymes)

    cutters   = {k: v for k, v in analysis.items() if v['cut_count'] > 0}
    nocutters = {k: v for k, v in analysis.items() if v['cut_count'] == 0}
    total_cuts = sum(v['cut_count'] for v in cutters.values())
    gcv = gc_content(seq_str)

    # ── Summary metrics ────────────────────────────────────────────────────
    st.markdown(f"""
    <div class="metric-row">
      <div class="metric-card">
        <div class="val">{len(seq_str):,}</div>
        <div class="lbl">Sequence length (bp)</div>
      </div>
      <div class="metric-card">
        <div class="val">{gcv}%</div>
        <div class="lbl">GC content</div>
      </div>
      <div class="metric-card">
        <div class="val">{len(cutters)}/{len(selected_enzymes)}</div>
        <div class="lbl">Enzymes that cut</div>
      </div>
      <div class="metric-card">
        <div class="val">{total_cuts}</div>
        <div class="lbl">Total cut sites</div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    # ── Enzyme summary badges ──────────────────────────────────────────────
    badge_html = '<div style="margin: 8px 0 16px;">'
    for name in selected_enzymes:
        data = analysis.get(name, {})
        if data.get('cuts_sequence'):
            badge_html += f'<span class="enzyme-badge cutter">{name} ✕{data["cut_count"]}</span>'
        else:
            badge_html += f'<span class="enzyme-badge nocutter">{name}</span>'
    badge_html += '</div>'
    st.markdown(badge_html, unsafe_allow_html=True)

    # ── Tabs ───────────────────────────────────────────────────────────────
    tabs = st.tabs(["CUT MAP", "VIRTUAL GEL", "FRAGMENTS", "GC PROFILE", "RAW DATA"])

    with tabs[0]:
        st.markdown('<p class="section-label">Restriction Map</p>', unsafe_allow_html=True)
        if cutters:
            plot_cut_map(len(seq_str), analysis, seq_str)
            st.markdown('<p class="section-label">Sequence with cut annotations</p>',
                        unsafe_allow_html=True)
            display_seq = format_sequence_display(seq_str[:600])
            st.markdown(f'<div class="seq-block">{display_seq}{"…" if len(seq_str)>600 else ""}</div>',
                        unsafe_allow_html=True)
        else:
            st.markdown('<div class="warn-box">None of the selected enzymes cut this sequence.</div>',
                        unsafe_allow_html=True)

    with tabs[1]:
        st.markdown('<p class="section-label">Simulated Gel Electrophoresis</p>',
                    unsafe_allow_html=True)
        st.markdown(
            '<div class="info-box">Fragment migration simulated on log scale. '
            'Left lane is a 100 bp DNA ladder. Hover bands for sizes.</div>',
            unsafe_allow_html=True
        )
        plot_gel(analysis, len(seq_str))

    with tabs[2]:
        st.markdown('<p class="section-label">Fragment Size Distribution</p>',
                    unsafe_allow_html=True)
        plot_fragment_distribution(analysis)

        st.markdown('<p class="section-label">Fragment Table</p>', unsafe_allow_html=True)
        rows = []
        for name, data in cutters.items():
            for i, frag in enumerate(data['fragments']):
                rows.append({
                    "Enzyme": name,
                    "Fragment": f"F{i+1}",
                    "Size (bp)": frag,
                    "Recognition Site": data['recognition_site'],
                    "Overhang": data['overhang'],
                })
        if rows:
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True,
                         column_config={
                             "Size (bp)": st.column_config.NumberColumn(format="%d bp"),
                         })
        else:
            st.info("No fragments to display.")

    with tabs[3]:
        st.markdown('<p class="section-label">GC Content Profile</p>', unsafe_allow_html=True)
        plot_gc_window(seq_str, gc_window)
        a_count = seq_str.upper().count('A')
        t_count = seq_str.upper().count('T')
        g_count = seq_str.upper().count('G')
        c_count = seq_str.upper().count('C')
        comp_df = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [a_count, t_count, g_count, c_count],
            "Percentage": [
                round(x / len(seq_str) * 100, 2)
                for x in [a_count, t_count, g_count, c_count]
            ]
        })
        fig_comp = go.Figure(go.Bar(
            x=comp_df["Nucleotide"], y=comp_df["Percentage"],
            marker=dict(
                color=['#4F7DF3', '#22C55E', '#14B8A6', '#F59E0B'],
                line=dict(width=0)
            ),
            hovertemplate="%{x}: %{y:.2f}%<extra></extra>"
        ))
        fig_comp.update_layout(
            **{k: v for k, v in PLOT_LAYOUT.items() if k not in ('xaxis', 'yaxis')},
            height=240,
            title=dict(text="Nucleotide Composition", font=dict(size=13, color='#F8FAFC'), x=0),
            xaxis=dict(title="Base", **PLOT_LAYOUT['xaxis']),
            yaxis=dict(title="%", **PLOT_LAYOUT['yaxis']),
        )
        st.plotly_chart(fig_comp, use_container_width=True)

    with tabs[4]:
        st.markdown('<p class="section-label">Enzyme Detail Table</p>', unsafe_allow_html=True)
        summary_rows = []
        for name in selected_enzymes:
            data = analysis.get(name, {})
            summary_rows.append({
                "Enzyme": name,
                "Cuts": data.get('cut_count', 0),
                "Recognition Site": data.get('recognition_site', '—'),
                "Overhang Type": data.get('overhang', '—'),
                "Cut Positions": ", ".join(str(p) for p in data.get('cuts', [])) or "—",
                "Fragment Sizes (bp)": ", ".join(str(f) for f in data.get('fragments', [])) or "—",
            })
        st.dataframe(pd.DataFrame(summary_rows), use_container_width=True, hide_index=True)

        st.markdown('<p class="section-label">Full Sequence</p>', unsafe_allow_html=True)
        full_display = format_sequence_display(seq_str)
        st.markdown(f'<div class="seq-block">{full_display}</div>', unsafe_allow_html=True)

        # FASTA export
        fasta_out = f">analyzed_sequence | {len(seq_str)} bp | GC={gcv}%\n"
        for i in range(0, len(seq_str), 60):
            fasta_out += seq_str[i:i+60] + "\n"
        st.download_button(
            label="Export FASTA",
            data=fasta_out,
            file_name="analyzed_sequence.fasta",
            mime="text/plain"
        )