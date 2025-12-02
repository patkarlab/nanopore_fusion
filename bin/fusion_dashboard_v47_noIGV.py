#!/usr/bin/env python3
"""
Enhanced Fusion Caller Dashboard Generator (v47)
Run with: python fusion_dashboard_v47_noIGV.py your_file.xlsx --cytoband cytoBand.txt --gtf Homo_sapiens.GRCh38.105.gtf.gz
"""

import pandas as pd
import json
import sys
import argparse
import gzip
from pathlib import Path
import urllib.request
import os
import io
from collections import defaultdict

parser = argparse.ArgumentParser(description='Generate interactive fusion caller dashboard')
parser.add_argument('input_file', help='Input Excel file')
parser.add_argument('-o', '--output', help='Output HTML', default=None)
parser.add_argument('--min-callers', type=int, default=1)
parser.add_argument('--min-reads', type=int, default=1)
parser.add_argument('--cytoband', help='cytoBand.txt file', default=None)
parser.add_argument('--gtf', help='GTF file (can be .gtf or .gtf.gz)', default=None)
args = parser.parse_args()

input_file = args.input_file
if not Path(input_file).exists():
    print(f"ERROR: File not found: {input_file}")
    sys.exit(1)

output_file = args.output if args.output else f'{Path(input_file).stem}_dashboard.html'
print(f"Input: {input_file}\nOutput: {output_file}\n")

# Load cytoband
cytoband_data = None
if args.cytoband and Path(args.cytoband).exists():
    try:
        cytoband_df = pd.read_csv(args.cytoband, sep='\t', names=['chrom','start','end','name','stain'])
        cytoband_data = cytoband_df.to_dict('records')
        print(f"Loaded {len(cytoband_data)} cytobands\n")
    except Exception as e:
        print(f"Warning: Could not load cytoband file: {e}")


def parse_gtf_file(gtf_path):
    """Parse GTF file and extract exon, transcript biotype, and transcript support level information."""
    gene_structures = {}

    if not gtf_path or not Path(gtf_path).exists():
        print("Warning: GTF file not found  using placeholder data.\n")
        return gene_structures

    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'exon':
                continue

            chrom, start, end, strand, attrs = (
                parts[0],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                parts[8],
            )

            gene_name = transcript_id = biotype = tsl = None
            exon_number = None

            for attr in attrs.split(';'):
                attr = attr.strip()
                if not attr:
                    continue

                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                elif attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                elif attr.startswith('exon_number'):
                    val = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    try:
                        exon_number = int(val)
                    except ValueError:
                        exon_number = None
                elif attr.startswith('transcript_support_level'):
                    val = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    # Normalize 1 (assigned to previous version &)
                    tsl = val.strip().split()[0].strip('()"')
                elif attr.startswith('transcript_biotype') or attr.startswith('transcript_type'):
                    biotype = attr.split('"')[1] if '"' in attr else attr.split()[-1]

            if not (gene_name and transcript_id and exon_number):
                continue

            g = gene_structures.setdefault(
                gene_name,
                {'chr': chrom, 'strand': strand, 'start': start, 'end': end, 'transcripts': {}},
            )
            g['start'] = min(g['start'], start)
            g['end'] = max(g['end'], end)

            t = g['transcripts'].setdefault(
                transcript_id, {'tsl': tsl, 'biotype': biotype, 'exons': []}
            )
            t['tsl'] = tsl or t.get('tsl')
            t['biotype'] = biotype or t.get('biotype')
            t['exons'].append({'exon_number': exon_number, 'start': start, 'end': end})

    # Sort exons within each transcript (respect strand)
    for g in gene_structures.values():
        for t in g['transcripts'].values():
            t['exons'].sort(key=lambda e: e['start'], reverse=(g['strand'] == '-'))

    return gene_structures


def tsl_val(t):
    """Return normalized TSL string (e.g. '1', '2', etc.)"""
    raw_tsl = t.get('tsl') or ''
    if not raw_tsl.strip():
        return None
    # Extract the first numeric or roman part (handles text like '1 (assigned...)')
    first_token = raw_tsl.strip().split()[0]
    return first_token.strip('()"')

def is_pc(t):
    """Check if transcript is protein_coding"""
    return (t.get('biotype') or '').lower() == 'protein_coding'

def find_best_transcript(gene_structure):
    """Select the best transcript: prefer protein_coding + TSL=1 (or '1 (assigned...)'), else longest"""
    transcripts = gene_structure.get('transcripts', {})
    if not transcripts:
        return None

    # Step 1: try to find TSL=1 and protein_coding
    pc_tsl1 = [
        tid for tid, t in transcripts.items()
        if is_pc(t) and tsl_val(t) and tsl_val(t).startswith('1')
    ]
    if pc_tsl1:
        # If multiple, choose the one with the most exons
        return max(pc_tsl1, key=lambda tid: len(transcripts[tid].get('exons', [])))

    # Step 2: fall back to protein_coding with max exons
    pc_transcripts = [tid for tid, t in transcripts.items() if is_pc(t)]
    if pc_transcripts:
        return max(pc_transcripts, key=lambda tid: len(transcripts[tid].get('exons', [])))

    # Step 3: fallback to transcript with most exons
    return max(transcripts, key=lambda tid: len(transcripts[tid].get('exons', [])))

def find_breakpoint_exon(gene_structure, position):
    """Find which exon contains or is closest to the breakpoint position."""
    if not gene_structure or 'transcripts' not in gene_structure:
        return None, None

    best_tid = find_best_transcript(gene_structure)
    if not best_tid:
        return None, None

    transcript_data = gene_structure['transcripts'][best_tid]
    exons = transcript_data.get('exons', [])
    if not exons:
        return best_tid, None

    for exon in exons:
        if exon['start'] <= position <= exon['end']:
            return best_tid, exon['exon_number']

    # If not inside any exon, find closest exon
    closest = min(exons, key=lambda e: min(abs(e['start'] - position), abs(e['end'] - position)))
    return best_tid, closest['exon_number']

gtf_structures = parse_gtf_file(args.gtf)

fusion_callers = {
    'jaffal': {'fusion_col':'fusion genes','reads_col':'spanning reads','display_name':'Jaffal',
               'chr1_col':'chrom1','pos1_col':'base1','chr2_col':'chrom2','pos2_col':'base2',
               'confidence_col':'classification','known_col':'known'},
    'longgf_final': {'fusion_col':'gene','reads_col':'num_reads','display_name':'LongGF',
                     'breakpoint1_col':'breakpoint1','breakpoint2_col':'breakpoint2'},
    'genion_final': {'fusion_col':'gene','reads_col':'num_reads','display_name':'Genion','ranges_col':'ranges'},
    'ctat-LR.fusion_predictions': {'fusion_col':'#FusionName','reads_col':'num_LR','display_name':'CTAT-LR',
                                    'chr1_col':'LeftBreakpoint','chr2_col':'RightBreakpoint','annots_col':'annots'},
    'fusionseeker_final': {'fusion_col':None,'reads_col':'NumSupp','display_name':'FusionSeeker',
                           'chr1_col':'Chrom1','pos1_col':'Breakpoint1','chr2_col':'Chrom2','pos2_col':'Breakpoint2'}
}

all_fusions = {}

def parse_bp(s):
    """Parse breakpoint string like 'chr1:12345' into chromosome and position"""
    if pd.isna(s) or s=='': 
        return None, None
    p = str(s).split(':')
    if len(p) >= 2:
        try: 
            return p[0], int(p[1])
        except: 
            return None, None
    return None, None

print("Processing Excel file...")
try:
    xls = pd.ExcelFile(input_file)
    sheets = xls.sheet_names
    print(f"Found {len(sheets)} sheets\n")
except Exception as e:
    print(f"ERROR reading Excel file: {e}")
    sys.exit(1)

# Map expected sheet names to actual sheet names
mapping = {}
for exp in fusion_callers:
    for act in sheets:
        # Special handling for ctat-lr variations
        if exp == 'ctat-LR.fusion_predictions':
            if 'ctat' in act.lower() and 'lr' in act.lower():
                mapping[exp] = act
                break
        elif exp.lower() in act.lower():
            mapping[exp] = act
            break

# Process each caller's data
for caller, cols in fusion_callers.items():
    if caller not in mapping: 
        print(f"  Skipping {caller} - sheet not found")
        continue
    
    sheet = mapping[caller]
    print(f"  Processing {sheet}...")
    
    try:
        df = pd.read_excel(input_file, sheet_name=sheet)
        if df.empty: 
            print(f"    Sheet is empty")
            continue
        
        print(f"    Found {len(df)} fusions")
        
        for _, row in df.iterrows():
            try:
                # Get fusion name
                if caller == 'fusionseeker_final':
                    fname = f"{row['Gene1']}:{row['Gene2']}"
                else:
                    fname = row[cols['fusion_col']]
                
                fname = str(fname).replace('::', ':').replace('--', ':')
                reads = int(row[cols['reads_col']])
                
                # Extract breakpoint information
                c1, p1, c2, p2 = None, None, None, None
                
                if caller == 'jaffal':
                    c1 = row[cols['chr1_col']]
                    p1 = int(row[cols['pos1_col']])
                    c2 = row[cols['chr2_col']]
                    p2 = int(row[cols['pos2_col']])
                    
                elif caller == 'longgf_final':
                    c1, p1 = parse_bp(row[cols['breakpoint1_col']])
                    c2, p2 = parse_bp(row[cols['breakpoint2_col']])
                    
                elif caller == 'genion_final':
                    rs = str(row.get(cols['ranges_col'], ''))
                    if rs != 'nan':
                        rng = rs.split(';')
                        if len(rng) >= 2:
                            pt1 = rng[0].split(':')
                            if len(pt1) >= 2: 
                                c1 = pt1[0]
                                p1 = int(pt1[1].split('-')[0])
                            pt2 = rng[1].split(':')
                            if len(pt2) >= 2: 
                                c2 = pt2[0]
                                p2 = int(pt2[1].split('-')[0])
                                
                elif caller == 'ctat-LR.fusion_predictions':
                    c1, p1 = parse_bp(row[cols['chr1_col']])
                    c2, p2 = parse_bp(row[cols['chr2_col']])
                    
                elif caller == 'fusionseeker_final':
                    c1 = row[cols['chr1_col']]
                    p1 = int(row[cols['pos1_col']])
                    c2 = row[cols['chr2_col']]
                    p2 = int(row[cols['pos2_col']])
                
                # Store fusion data
                if fname not in all_fusions:
                    all_fusions[fname] = {
                        'fusion': fname,
                        'callers': {},
                        'breakpoints': {},
                        'caller_max_reads': {},  # Track max reads per caller
                        'chr1': c1,
                        'pos1': p1,
                        'chr2': c2,
                        'pos2': p2,
                        'confidence': None,
                        'known': None,
                        'annotations': '[]',
                        'max_jaffal_reads': 0
                    }
                
                # Only update reads and breakpoints if this entry has more reads than previous entries for this caller
                caller_name = cols['display_name']
                if caller_name not in all_fusions[fname]['caller_max_reads'] or reads > all_fusions[fname]['caller_max_reads'][caller_name]:
                    all_fusions[fname]['caller_max_reads'][caller_name] = reads
                    all_fusions[fname]['callers'][caller_name] = reads
                    
                    if c1 and p1 and c2 and p2:
                        all_fusions[fname]['breakpoints'][caller_name] = {
                            'chr1': c1, 'pos1': p1, 'chr2': c2, 'pos2': p2
                        }
                
                # Store additional metadata - for jaffal, only use data from the entry with highest spanning reads
                if caller == 'jaffal':
                    if reads > all_fusions[fname]['max_jaffal_reads']:
                        all_fusions[fname]['max_jaffal_reads'] = reads
                        all_fusions[fname]['confidence'] = row.get(cols.get('confidence_col', ''), 'Unknown')
                        all_fusions[fname]['known'] = 'Yes' if row.get(cols.get('known_col', ''), '-') == 'Yes' else 'No'
                
                if caller == 'ctat-LR.fusion_predictions' and 'annots_col' in cols:
                    all_fusions[fname]['annotations'] = str(row.get(cols['annots_col'], '[]'))
                    
            except Exception as e:
                # Print parse error for visibility (was silently continuing)
                # Comment out the print if it becomes too verbose
                print(f"    Row parse error ({sheet}): {e}")
                continue
                
    except Exception as e:
        print(f"    Error processing sheet: {e}")
        continue

# Ensure all callers are represented
for fn in all_fusions:
    for c in ['Jaffal', 'LongGF', 'Genion', 'CTAT-LR', 'FusionSeeker']:
        if c not in all_fusions[fn]['callers']:
            all_fusions[fn]['callers'][c] = None

fusion_list = list(all_fusions.values())

#-----------------------------------------------
# Derive sample prefix from input filename (to be used for log file name)
#-----------------------------------------------
input_name = Path(args.input_file).stem  # removes directory and extension
log_file_path = f"{input_name}_gene_transcript_mapping.txt"

#----------------------------------------------
# === Utility: get all valid transcripts for dropdown
#----------------------------------------------
def get_all_valid_transcripts(g_struct):
    """Return all protein_coding + TSL1 transcripts, or fallback to all."""
    valid = []
    for tid, tdata in g_struct['transcripts'].items():
        biotype = (tdata.get('biotype') or '').lower()
        tsl_raw = (tdata.get('tsl') or '').strip()
        tsl_val = tsl_raw.split()[0].strip('()"') if tsl_raw else ''

        # accept protein_coding + TSL=1 (including "1 (assigned...)" cases)
        if 'protein_coding' in biotype and (tsl_val == '1' or tsl_val.startswith('1')):
            valid.append(tid)

    # fallback if no "valid" ones found
    if not valid:
        valid = list(g_struct['transcripts'].keys())
    return valid


# === Build gene structure data for each fusion ===
gene_structure_data = {}
if gtf_structures:
    print("Preparing gene structure data for fusions...\n")

    # derive prefix for transcript mapping file name
    input_prefix = Path(args.input_file).stem 
    mapping_out = f"{input_prefix}_gene_transcript_mapping.txt"

    with open(mapping_out, "w") as mapfile:
        mapfile.write("Fusion\tGene\tTranscript_ID\tBiotype\tTSL\tTotal_Exons\n")

        for fusion in fusion_list:
            genes = fusion['fusion'].split(':')
            if len(genes) != 2:
                continue

            gene1, gene2 = genes
            gene_structure_data[fusion['fusion']] = {'gene1': None, 'gene2': None}

            for side, gene_name, pos_key in [
                ('gene1', gene1, 'pos1'),
                ('gene2', gene2, 'pos2'),
            ]:
                if gene_name not in gtf_structures:
                    continue

                g_struct = gtf_structures[gene_name]
                transcript_id, breakpoint_exon = find_breakpoint_exon(g_struct, fusion.get(pos_key))
                transcript_id = transcript_id or 'N/A'

                available_tx = get_all_valid_transcripts(g_struct)

                # get details for the chosen transcript
                tx = g_struct['transcripts'].get(transcript_id, {})
                biotype = tx.get('biotype', 'N/A')
                tsl = tx.get('tsl', 'N/A')
                exons = [e['exon_number'] for e in tx.get('exons', [])]

                # write to mapping file
                mapfile.write(f"{fusion['fusion']}\t{gene_name}\t{transcript_id}\t{biotype}\t{tsl}\t{len(exons)}\n")

                # Store gene structure info
                # Store exon data for *all* transcripts of this gene
                all_tx_data = {}
                for tid, tdata in g_struct['transcripts'].items():
                    all_tx_data[tid] = {
                        'exons': tdata.get('exons', []),
                        'biotype': tdata.get('biotype', 'N/A'),
                        'tsl': tdata.get('tsl', 'N/A')
                    }
                
                # Store gene structure info
                gene_structure_data[fusion['fusion']][side] = {
                    'name': gene_name,
                    'chr': g_struct['chr'],
                    'pos': fusion.get(pos_key),
                    'strand': g_struct['strand'],
                    'exons': exons,
                    'total_exons': len(exons),
                    'breakpoint_exon': breakpoint_exon,
                    'transcript': transcript_id,
                    'available_transcripts': available_tx,
                    # include exon info for all transcripts
                    'all_transcripts': {
                        tid: {
                            'exons': t.get('exons', []),
                            'biotype': t.get('biotype', 'N/A'),
                            'tsl': t.get('tsl', 'N/A'),
                        }
                        for tid, t in g_struct['transcripts'].items()
                    }
                }



    print(f"\n Prepared structure data for {len(gene_structure_data)} fusions.")
    print(f" Transcript mapping written to: {mapping_out}\n")



# Load coverage data if available
coverage_data = None
bedtools_data = None

if any('mosdepth' in s.lower() for s in sheets):
    try:
        ms = next((s for s in sheets if 'mosdepth' in s.lower()), None)
        if ms:
            print(f"\nLoading coverage (mosdepth): {ms}")
            coverage_data = pd.read_excel(input_file, sheet_name=ms).to_dict('records')
            print(f"  Loaded {len(coverage_data)} regions")
    except Exception as e:
        print(f"  Warning: Could not load mosdepth coverage data: {e}")

if any('bedtools' in s.lower() for s in sheets):
    try:
        bt = next((s for s in sheets if 'bedtools' in s.lower()), None)
        if bt:
            print(f"Loading coverage (bedtools): {bt}")
            bedtools_data = pd.read_excel(input_file, sheet_name=bt).to_dict('records')
            print(f"  Loaded {len(bedtools_data)} regions")
    except Exception as e:
        print(f"  Warning: Could not load bedtools coverage data: {e}")

print(f"\n Total fusions found: {len(fusion_list)}")

if len(fusion_list) == 0:
    print("ERROR: No fusions found in the input file!")
    # Instead of exiting immediately, produce a minimal HTML indicating no data
    # (This mirrors behavior requested: avoid blank output)
    html_empty = f"""<!DOCTYPE html>
    <html><head><title>Fusion Dashboard - {Path(input_file).stem}</title></head><body>
    <h1>No fusions found in {Path(input_file).name}</h1>
    <p>Check the sheet names and column headers. Script run completed.</p>
    </body></html>"""
    with open(output_file, 'w') as fh:
        fh.write(html_empty)
    print(f"Wrote minimal HTML to {output_file}")
    sys.exit(1)

print(f"\nGenerating HTML dashboard...")

# Generate HTML
sample_name = Path(input_file).stem
sel3 = 'selected' if args.min_callers == 3 else ''
sel4 = 'selected' if args.min_callers == 4 else ''
sel5 = 'selected' if args.min_callers == 5 else ''

# Create the complete HTML file
html_content = f"""<!DOCTYPE html>
<html><head><title>Fusion Dashboard - {sample_name}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); min-height: 100vh; padding: 20px; }}
.dashboard {{ max-width: 1600px; margin: 0 auto; }}
.header {{ background: white; padding: 25px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 20px; }}
.header h1 {{ color: #2c3e50; font-size: 32px; margin-bottom: 10px; }}
.header p {{ color: #7f8c8d; font-size: 16px; }}
.controls {{ background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 20px; display: flex; gap: 20px; flex-wrap: wrap; }}
.control-group {{ display: flex; flex-direction: column; gap: 5px; }}
.control-group label {{ font-weight: 600; color: #2c3e50; font-size: 14px; }}
.control-group select, .control-group input {{ padding: 8px 12px; border: 2px solid #e0e0e0; border-radius: 6px; font-size: 14px; min-width: 150px; }}
.stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 20px; }}
.stat-card {{ background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
.stat-card h3 {{ color: #7f8c8d; font-size: 14px; font-weight: 600; margin-bottom: 8px; }}
.stat-card .value {{ color: #2c3e50; font-size: 32px; font-weight: bold; }}
.main-content {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
.panel {{ background: white; padding: 25px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
.panel h2 {{ color: #2c3e50; font-size: 20px; margin-bottom: 15px; padding-bottom: 10px; border-bottom: 2px solid #e0e0e0; }}
.table-container {{ overflow-x: auto; margin-top: 15px; max-height: 600px; overflow-y: auto; }}
table {{ width: 100%; border-collapse: collapse; }}
th {{ background-color: #667eea; color: white; padding: 12px; text-align: left; font-weight: 600; position: sticky; top: 0; z-index: 10; }}
td {{ padding: 10px 12px; border-bottom: 1px solid #e0e0e0; }}
tr:hover {{ background-color: #f8f9fa; cursor: pointer; }}
.fusion-name {{ font-weight: 600; color: #2c3e50; }}
.read-count {{ text-align: center; font-weight: 600; color: #27ae60; }}
.no-data {{ text-align: center; color: #bdc3c7; }}
#circos {{ width: 100%; height: 700px; display: flex; align-items: center; justify-content: center; }}
.fusion-arc {{ fill: none; stroke-width: 2; opacity: 0.6; cursor: pointer; transition: all 0.3s; }}
.fusion-arc:hover {{ opacity: 1; stroke-width: 3; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.9); color: white; padding: 12px; border-radius: 6px; font-size: 13px; pointer-events: none; opacity: 0; z-index: 1000; transition: opacity 0.2s; }}
.chr-label {{ font-size: 11px; font-weight: bold; fill: #2c3e50; }}
.breakpoint-marker {{ fill: #e74c3c; stroke: white; stroke-width: 2; cursor: pointer; }}
.details-panel {{ background: white; padding: 25px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-top: 20px; display: none; }}
.details-panel.active {{ display: block; }}
.detail-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin-top: 15px; }}
.detail-item {{ padding: 15px; background: #f8f9fa; border-radius: 8px; border-left: 4px solid #667eea; }}
.detail-item h4 {{ color: #7f8c8d; font-size: 12px; font-weight: 600; margin-bottom: 5px; text-transform: uppercase; }}
.detail-item p {{ color: #2c3e50; font-size: 16px; font-weight: 600; }}
.breakpoint-table {{ width: 100%; margin-top: 15px; }}
.breakpoint-table th {{ background-color: #667eea; color: white; padding: 10px; text-align: left; }}
.breakpoint-table td {{ padding: 8px 10px; border-bottom: 1px solid #e0e0e0; }}
.close-details {{ float: right; background: #e74c3c; color: white; border: none; padding: 8px 15px; border-radius: 6px; cursor: pointer; font-weight: 600; }}
.close-details:hover {{ background: #c0392b; }}
.close-details:hover {{ background: #c0392b; }}
/* === Transcript selector styles === */
#transcriptSelectors {{
  margin-top: 10px;
  margin-bottom: 15px;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 10px;
}}
#transcriptSelectors label {{
  font-weight: 600;
  color: #2c3e50;
}}
#transcriptSelectors select {{
  background: #f8f9fa;
  padding: 6px 10px;
  border-radius: 6px;
  border: 1.5px solid #ccc;
  font-size: 13px;
  transition: border-color 0.2s ease;
}}
#transcriptSelectors select:focus {{
  outline: none;
  border-color: #667eea;
}}
.coverage-tabs, .breakpoint-tabs {{ display: flex; gap: 10px; margin-bottom: 15px; }}
.coverage-tab, .breakpoint-tab {{ 
  padding: 8px 20px; 
  background: #f8f9fa; 
  border: 2px solid #e0e0e0; 
  border-radius: 6px; 
  cursor: pointer; 
  font-weight: 600; 
  color: #7f8c8d;
  transition: all 0.3s;
}}
.coverage-tab:hover, .breakpoint-tab:hover {{ background: #e0e0e0; }}
.coverage-tab.active, .breakpoint-tab.active {{ 
  background: #667eea; 
  color: white; 
  border-color: #667eea;
}}
.read-alignment {{ 
  font-family: 'Courier New', monospace; 
  font-size: 11px; 
  padding: 4px 8px; 
  margin: 2px 0;
  background: #f8f9fa;
  border-left: 3px solid #3498db;
  white-space: nowrap;
  overflow-x: auto;
}}
.read-alignment.reverse {{ border-left-color: #e74c3c; }}
.breakpoint-position {{ background: #ffeb3b; font-weight: bold; }}
</style></head><body>
<div class="dashboard">
<div class="header"><h1>>Fusion Caller Interactive Dashboard</h1><p>Sample: {sample_name}</p></div>
<div class="controls">
<div class="control-group"><label>Minimum Callers</label><select id="minCallers">
<option value="1">1+</option><option value="2">2+</option><option value="3" {sel3}>3+</option>
<option value="4" {sel4}>4+</option><option value="5" {sel5}>5</option></select></div>
<div class="control-group"><label>Minimum Reads</label><input type="number" id="minReads" value="{args.min_reads}" min="0"></div>
<div class="control-group"><label>Sort By</label><select id="sortBy">
<option value="callers">Callers</option><option value="reads">Reads</option><option value="name">Name</option></select></div>
</div>
<div class="stats">
<div class="stat-card"><h3>Total Fusions</h3><div class="value" id="totalFusions">0</div></div>
<div class="stat-card"><h3>High Confidence (40+ reads)</h3><div class="value" id="highConfidence">0</div></div>
<div class="stat-card"><h3>Average Reads per fusion</h3><div class="value" id="avgReads">0</div></div>
<div class="stat-card"><h3>Chromosomes</h3><div class="value" id="chrInvolved">0</div></div>
</div>
<div class="main-content">
<div class="panel"><h2>Fusion Table</h2><div class="table-container">
<table id="fusionTable"><thead><tr><th>Fusion</th><th>Jaffal</th><th>LongGF</th><th>Genion</th>
<th>CTAT-LR</th><th>FusionSeeker</th><th>Callers</th></tr></thead><tbody id="tableBody"></tbody></table>
</div></div>
<div class="panel"><h2>Circos Plot</h2><div id="circos"></div></div>
</div>
<div class="details-panel" id="detailsPanel">
<h2><span id="detailsFusionName">Fusion Details</span><button class="close-details" onclick="closeDetails()">Close</button></h2>
<div class="detail-grid" id="detailsGrid"></div>
<h3 style="margin-top:20px;color:#2c3e50;">Breakpoint Information</h3>
<table class="breakpoint-table" id="breakpointTable"><thead><tr>
<th>Caller</th><th>Chromosome 1</th><th>Position 1</th><th>Chromosome 2</th><th>Position 2</th>
</tr></thead><tbody id="breakpointTableBody"></tbody></table>
    <div id="coverageSection" style="display:none;margin-top:20px;">
<h3 style="color:#2c3e50;">Gene Fusion Diagram</h3>
<div id="transcriptSelectors" style="margin:10px 0 15px 0;">
    <label>Gene 1 transcript:</label>
    <select id="gene1TranscriptSelect"
            onchange="updateFusionTranscript('gene1', this.value)"
            style="margin-left:6px; padding:4px 6px; border-radius:4px; border:1px solid #ccc; font-size:12px;"></select>
    <label style="margin-left:20px;">Gene 2 transcript:</label>
    <select id="gene2TranscriptSelect"
            onchange="updateFusionTranscript('gene2', this.value)"
            style="margin-left:6px; padding:4px 6px; border-radius:4px; border:1px solid #ccc; font-size:12px;"></select>
</div>
<div id="geneFusionDiagram" style="width:100%;height:400px;margin-bottom:30px;"></div>
<h3 style="color:#2c3e50;">Coverage Plot</h3>
<div id="coveragePlot" style="width:100%;height:400px;"></div>
<h3 style="color:#2c3e50;margin-top:20px;">Coverage Details</h3>
<div class="coverage-tabs" style="margin-top:10px;">
    <button class="coverage-tab active" onclick="switchCoverageTab('mosdepth', event)">Mosdepth</button>
    <button class="coverage-tab" onclick="switchCoverageTab('bedtools', event)">Bedtools</button>
</div>
<div id="mosdepthTable" class="coverage-table-container" style="width:100%;margin-top:10px;overflow-x:auto;display:block;"></div>
<div id="bedtoolsTable" class="coverage-table-container" style="width:100%;margin-top:10px;overflow-x:auto;display:none;"></div>
</div>
</div>
</div>
<div class="tooltip" id="tooltip"></div>
<script>
const fusionData = {json.dumps(fusion_list)};
const cytobandData = {json.dumps(cytoband_data if cytoband_data else [])};
const coverageData = {json.dumps(coverage_data if coverage_data else [])};
const bedtoolsData = {json.dumps(bedtools_data if bedtools_data else [])};
const geneStructureData = {json.dumps(gene_structure_data)};
const chromosomes = [
{{name:'chr1',length:248956422,color:'#e74c3c'}},{{name:'chr2',length:242193529,color:'#3498db'}},
{{name:'chr3',length:198295559,color:'#2ecc71'}},{{name:'chr4',length:190214555,color:'#f39c12'}},
{{name:'chr5',length:181538259,color:'#9b59b6'}},{{name:'chr6',length:170805979,color:'#1abc9c'}},
{{name:'chr7',length:159345973,color:'#e67e22'}},{{name:'chr8',length:145138636,color:'#34495e'}},
{{name:'chr9',length:138394717,color:'#16a085'}},{{name:'chr10',length:133797422,color:'#27ae60'}},
{{name:'chr11',length:135086622,color:'#2980b9'}},{{name:'chr12',length:133275309,color:'#8e44ad'}},
{{name:'chr13',length:114364328,color:'#c0392b'}},{{name:'chr14',length:107043718,color:'#d35400'}},
{{name:'chr15',length:101991189,color:'#2c3e50'}},{{name:'chr16',length:90338345,color:'#7f8c8d'}},
{{name:'chr17',length:83257441,color:'#c0392b'}},{{name:'chr18',length:80373285,color:'#16a085'}},
{{name:'chr19',length:58617616,color:'#d35400'}},{{name:'chr20',length:64444167,color:'#27ae60'}},
{{name:'chr21',length:46709983,color:'#8e44ad'}},{{name:'chr22',length:50818468,color:'#2980b9'}},
{{name:'chrX',length:156040895,color:'#e91e63'}},{{name:'chrY',length:57227415,color:'#9c27b0'}}
];

let filteredData = [...fusionData];
let selectedFusion = null;

function updateTable() {{
    const minC = parseInt(document.getElementById('minCallers').value);
    const minR = parseInt(document.getElementById('minReads').value);
    const sortBy = document.getElementById('sortBy').value;
    
    filteredData = fusionData.filter(f => {{
        const callerCount = Object.values(f.callers).filter(v => v !== null).length;
        const totalReads = Object.values(f.callers).filter(v => v !== null).reduce((a, b) => a + b, 0);
        return callerCount >= minC && totalReads >= minR;
    }});
    
    filteredData.sort((a, b) => {{
        if (sortBy === 'callers') {{
            return Object.values(b.callers).filter(v => v !== null).length - 
                   Object.values(a.callers).filter(v => v !== null).length;
        }}
        if (sortBy === 'reads') {{
            const aReads = Object.values(a.callers).filter(v => v !== null).reduce((sum, val) => sum + val, 0);
            const bReads = Object.values(b.callers).filter(v => v !== null).reduce((sum, val) => sum + val, 0);
            return bReads - aReads;
        }}
        return a.fusion.localeCompare(b.fusion);
    }});
    
    const tbody = document.getElementById('tableBody');
    tbody.innerHTML = '';
    
    filteredData.forEach(f => {{
        const row = document.createElement('tr');
        const callerCount = Object.values(f.callers).filter(v => v !== null).length;
        row.onclick = () => showFusionDetails(f);
        if (selectedFusion === f.fusion) {{
            row.style.backgroundColor = '#e3f2fd';
        }}
        row.innerHTML = `
            <td class="fusion-name">${{f.fusion}}</td>
            <td class="${{f.callers.Jaffal ? 'read-count' : 'no-data'}}">${{f.callers.Jaffal || '-'}}</td>
            <td class="${{f.callers.LongGF ? 'read-count' : 'no-data'}}">${{f.callers.LongGF || '-'}}</td>
            <td class="${{f.callers.Genion ? 'read-count' : 'no-data'}}">${{f.callers.Genion || '-'}}</td>
            <td class="${{f.callers['CTAT-LR'] ? 'read-count' : 'no-data'}}">${{f.callers['CTAT-LR'] || '-'}}</td>
            <td class="${{f.callers.FusionSeeker ? 'read-count' : 'no-data'}}">${{f.callers.FusionSeeker || '-'}}</td>
            <td class="read-count">${{callerCount}}</td>
        `;
        tbody.appendChild(row);
    }});
    
    document.getElementById('totalFusions').textContent = filteredData.length;
    
    const highConf = fusionData.filter(f => 
        Object.values(f.callers).filter(v => v !== null).reduce((a, b) => a + b, 0) >= 40
    ).length;
    document.getElementById('highConfidence').textContent = highConf;
    
    const avgReads = filteredData.length > 0 ? 
        Math.round(filteredData.reduce((sum, f) => 
            sum + Object.values(f.callers).filter(v => v !== null).reduce((a, b) => a + b, 0), 0
        ) / filteredData.length) : 0;
    document.getElementById('avgReads').textContent = avgReads;
    
    const chrs = new Set();
    filteredData.forEach(f => {{
        if (f.chr1) chrs.add(f.chr1);
        if (f.chr2) chrs.add(f.chr2);
    }});
    document.getElementById('chrInvolved').textContent = chrs.size;
    
    drawCircos();
}}

function drawCircos() {{
    const svg = d3.select('#circos');
    svg.selectAll('*').remove();
    
    const width = 700;
    const height = 700;
    const radius = Math.min(width, height) / 2 - 100;
    
    const g = svg.append('svg')
        .attr('width', width)
        .attr('height', height)
        .append('g')
        .attr('transform', `translate(${{width/2}},${{height/2}})`);
    
    const totalLength = chromosomes.reduce((sum, chr) => sum + chr.length, 0);
    let currentAngle = 0;
    const chrData = [];
    
    chromosomes.forEach(chr => {{
        const angleSpan = (chr.length / totalLength) * 2 * Math.PI;
        chrData.push({{
            ...chr,
            startAngle: currentAngle,
            endAngle: currentAngle + angleSpan
        }});
        currentAngle += angleSpan;
    }});
    
    // Draw cytoband patterns if available
    if (cytobandData.length > 0) {{
        const bandArc = d3.arc()
            .innerRadius(radius - 20)
            .outerRadius(radius);
        
        cytobandData.forEach(band => {{
            const chr = chrData.find(c => c.name === band.chrom);
            if (!chr) return;
            
            const bandStart = chr.startAngle + (band.start / chr.length) * (chr.endAngle - chr.startAngle);
            const bandEnd = chr.startAngle + (band.end / chr.length) * (chr.endAngle - chr.startAngle);
            
            let bandColor = chr.color;
            let bandOpacity = 0.3;
            
            if (band.stain === 'gneg') {{ bandOpacity = 0.2; }}
            else if (band.stain === 'gpos25') {{ bandOpacity = 0.4; }}
            else if (band.stain === 'gpos50') {{ bandOpacity = 0.6; }}
            else if (band.stain === 'gpos75') {{ bandOpacity = 0.8; }}
            else if (band.stain === 'gpos100') {{ bandOpacity = 1.0; }}
            else if (band.stain === 'acen') {{ bandColor = '#d62728'; bandOpacity = 0.7; }}
            else if (band.stain === 'gvar') {{ bandOpacity = 0.5; }}
            else if (band.stain === 'stalk') {{ bandOpacity = 0.3; }}
            
            g.append('path')
                .attr('d', bandArc({{startAngle: bandStart, endAngle: bandEnd}}))
                .attr('fill', bandColor)
                .attr('opacity', bandOpacity)
                .attr('stroke', 'white')
                .attr('stroke-width', 0.5);
        }});
    }} else {{
        const arc = d3.arc()
            .innerRadius(radius - 20)
            .outerRadius(radius);
        
        g.selectAll('.chr-arc')
            .data(chrData)
            .enter()
            .append('path')
            .attr('class', 'chr-arc')
            .attr('d', arc)
            .attr('fill', d => d.color)
            .attr('opacity', 0.7);
    }}
    
    // Chromosome labels with collision avoidance
    const labelRadius = radius + 25;
    const labels = [];
    
    chrData.forEach(chr => {{
        const angle = (chr.startAngle + chr.endAngle) / 2;
        const x = labelRadius * Math.cos(angle - Math.PI / 2);
        const y = labelRadius * Math.sin(angle - Math.PI / 2);
        labels.push({{chr: chr.name.replace('chr', ''), x, y, angle}});
    }});
    
    g.selectAll('.chr-label')
        .data(labels)
        .enter()
        .append('text')
        .attr('class', 'chr-label')
        .attr('transform', d => `translate(${{d.x}},${{d.y}})`)
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'middle')
        .style('font-size', '11px')
        .text(d => d.chr);
    
    const tooltip = d3.select('#tooltip');
    
    // Sort fusions by total reads for better layering
    const sortedFusions = [...filteredData].sort((a, b) => {{
        const aReads = Object.values(a.callers).filter(v => v !== null).reduce((sum, val) => sum + val, 0);
        const bReads = Object.values(b.callers).filter(v => v !== null).reduce((sum, val) => sum + val, 0);
        return aReads - bReads;
    }});
    
    sortedFusions.forEach(fusion => {{
        if (!fusion.chr1 || !fusion.chr2 || !fusion.pos1 || !fusion.pos2) return;
        
        const chr1Data = chrData.find(c => c.name === fusion.chr1);
        const chr2Data = chrData.find(c => c.name === fusion.chr2);
        
        if (!chr1Data || !chr2Data) return;
        
        const angle1 = chr1Data.startAngle + (fusion.pos1 / chr1Data.length) * (chr1Data.endAngle - chr1Data.startAngle);
        const angle2 = chr2Data.startAngle + (fusion.pos2 / chr2Data.length) * (chr2Data.endAngle - chr2Data.startAngle);
        
        const x1 = (radius - 30) * Math.cos(angle1 - Math.PI / 2);
        const y1 = (radius - 30) * Math.sin(angle1 - Math.PI / 2);
        const x2 = (radius - 30) * Math.cos(angle2 - Math.PI / 2);
        const y2 = (radius - 30) * Math.sin(angle2 - Math.PI / 2);
        
        const totalReads = Object.values(fusion.callers).filter(v => v !== null).reduce((a, b) => a + b, 0);
        const strokeWidth = Math.max(1.5, Math.min(5, totalReads / 20));
        
        const arcPath = g.append('path')
            .attr('class', 'fusion-arc')
            .attr('data-fusion', fusion.fusion)
            .attr('d', `M ${{x1}} ${{y1}} Q 0 0 ${{x2}} ${{y2}}`)
            .attr('stroke', chr1Data.color)
            .attr('stroke-width', strokeWidth)
            .attr('opacity', selectedFusion === fusion.fusion ? 1 : 0.6)
            .on('mouseover', function(event) {{
                if (selectedFusion !== fusion.fusion) {{
                    d3.select(this).attr('opacity', 0.9).attr('stroke-width', strokeWidth + 2);
                }}
                tooltip.style('opacity', 1)
                    .html(`<strong>${{fusion.fusion}}</strong><br>Total Reads: ${{totalReads}}<br>Callers: ${{Object.values(fusion.callers).filter(v => v !== null).length}}`)
                    .style('left', (event.pageX + 10) + 'px')
                    .style('top', (event.pageY - 28) + 'px');
            }})
            .on('mouseout', function() {{
                if (selectedFusion !== fusion.fusion) {{
                    d3.select(this).attr('opacity', 0.6).attr('stroke-width', strokeWidth);
                }}
                tooltip.style('opacity', 0);
            }})
            .on('click', () => showFusionDetails(fusion));
        
        // Draw breakpoint markers
        g.append('circle')
            .attr('class', 'breakpoint-marker')
            .attr('cx', x1)
            .attr('cy', y1)
            .attr('r', 4)
            .style('cursor', 'pointer')
            .on('click', () => showFusionDetails(fusion));
        
        g.append('circle')
            .attr('class', 'breakpoint-marker')
            .attr('cx', x2)
            .attr('cy', y2)
            .attr('r', 4)
            .style('cursor', 'pointer')
            .on('click', () => showFusionDetails(fusion));
    }});
}}

function showFusionDetails(fusion) {{
    selectedFusion = fusion.fusion;
    
    // Highlight selected fusion in circos
    d3.selectAll('.fusion-arc').attr('opacity', function() {{
        return d3.select(this).attr('data-fusion') === fusion.fusion ? 1 : 0.3;
    }}).attr('stroke-width', function() {{
        const currentWidth = parseFloat(d3.select(this).attr('stroke-width'));
        return d3.select(this).attr('data-fusion') === fusion.fusion ? currentWidth + 2 : currentWidth;
    }});
    
    // Highlight in table
    document.querySelectorAll('#tableBody tr').forEach(row => {{
        row.style.backgroundColor = '';
    }});
    const rows = document.querySelectorAll('#tableBody tr');
    rows.forEach(row => {{
        if (row.cells[0].textContent === fusion.fusion) {{
            row.style.backgroundColor = '#e3f2fd';
        }}
    }});
    
    document.getElementById('detailsFusionName').textContent = fusion.fusion;
    
    const detailsGrid = document.getElementById('detailsGrid');
    detailsGrid.innerHTML = '';
    
    const totalReads = Object.values(fusion.callers).filter(v => v !== null).reduce((a, b) => a + b, 0);
    const callerCount = Object.values(fusion.callers).filter(v => v !== null).length;
    
    // Parse annots data (ctat-lr 'annots' -> stored as 'annotations' in fusion JSON)
    let reportedIn = 'N/A';
    let fusionType = 'N/A';
    
    const rawAnnots = (fusion.annotations !== undefined && fusion.annotations !== null) ? fusion.annotations : fusion.annots;
    if (rawAnnots) {{
        try {{
            let arr = null;
            if (Array.isArray(rawAnnots)) {{
                arr = rawAnnots;
            }} else if (typeof rawAnnots === 'string') {{
                const s = rawAnnots.trim();
                if (s.startsWith('[') && s.endsWith(']')) {{
                    arr = JSON.parse(s);
                }} else {{
                    // Fallback: split CSV-like string
                    const inner = s.replace(/^\"|\"$/g, '');
                    arr = inner.split(',').map(x => x.trim().replace(/^\"|\"$/g, '').replace(/^'|'$/g, ''));
                }}
            }}
            if (arr && Array.isArray(arr) && arr.length > 0) {{
                const last = String(arr[arr.length - 1]).trim();
                fusionType = last || 'N/A';
                if (arr.length > 1) {{
                    const dbs = arr.slice(0, -1).map(x => String(x).trim()).filter(Boolean);
                    if (dbs.length > 0) reportedIn = dbs.join(', ');
                }}
            }}
        }} catch (e) {{
            console.log('Could not parse annotations:', e, rawAnnots);
        }}
    }}
    
    detailsGrid.innerHTML = `
        <div class="detail-item">
            <h4>Classification</h4>
            <p>${{fusion.confidence || 'N/A'}}</p>
        </div>
        <div class="detail-item">
            <h4>Known Fusion</h4>
            <p>${{fusion.known || 'N/A'}}</p>
        </div>
        <div class="detail-item">
            <h4>Reported In</h4>
            <p style="font-size: 11px; line-height: 1.4;">${{reportedIn}}</p>
        </div>
        <div class="detail-item">
            <h4>Callers Detected</h4>
            <p>${{callerCount}} / 5</p>
        </div>
        <div class="detail-item">
            <h4>Total Reads</h4>
            <p>${{totalReads}}</p>
        </div>
        <div class="detail-item">
            <h4>Fusion Type</h4>
            <p>${{fusionType}}</p>
        </div>
    `;
    
    const breakpointTableBody = document.getElementById('breakpointTableBody');
    breakpointTableBody.innerHTML = '';
    
    Object.entries(fusion.breakpoints).forEach(([caller, bp]) => {{
        const row = document.createElement('tr');
        row.innerHTML = `
            <td>${{caller}}</td>
            <td>${{bp.chr1}}</td>
            <td>${{bp.pos1.toLocaleString()}}</td>
            <td>${{bp.chr2}}</td>
            <td>${{bp.pos2.toLocaleString()}}</td>
        `;
        breakpointTableBody.appendChild(row);
    }});
    
    // Always show coverage section and tables
    const genes = fusion.fusion.split(':');
    const gene1 = genes[0];
    const gene2 = genes[1];
    
    document.getElementById('coverageSection').style.display = 'block';
    
    populateTranscriptDropdowns(fusion);
    // Draw fusion diagram
    drawGeneFusionDiagram(fusion, gene1, gene2);
    
    // Draw coverage plot if data available
    if (coverageData && coverageData.length > 0) {{
        const relevantCoverage = coverageData.filter(cov => {{
            const region = cov.region || '';
            return region.includes(gene1) || region.includes(gene2);
        }});
        
        if (relevantCoverage.length > 0) {{
            drawCoveragePlot(relevantCoverage, gene1, gene2);
        }} else {{
            document.getElementById('coveragePlot').innerHTML = '<p style="text-align:center;padding:40px;color:#7f8c8d;">No coverage data available for this fusion</p>';
        }}
    }} else {{
        document.getElementById('coveragePlot').innerHTML = '<p style="text-align:center;padding:40px;color:#7f8c8d;">No coverage data available</p>';
    }}

    
    // Always show coverage tables with all data
    drawMosdepthTable(coverageData);
    drawBedtoolsTable(bedtoolsData);
    
    document.getElementById('detailsPanel').classList.add('active');
    document.getElementById('detailsPanel').scrollIntoView({{ behavior: 'smooth' }});
}}

function drawGeneFusionDiagram(fusion, gene1, gene2) {{
    const diagramDiv = d3.select('#geneFusionDiagram');
    diagramDiv.selectAll('*').remove();

    // Get gene structure data from embedded data
    const geneData = geneStructureData[fusion.fusion] || {{ gene1: null, gene2: null }};

    // --- Responsive width and margins ---
    const containerWidth = document.getElementById('geneFusionDiagram').clientWidth || 1200;

    // Reduce left margin if screen is wide, keep some padding otherwise
    const leftMargin = containerWidth > 1000 ? 20 : 5;
    const rightMargin = 60;

    const width = containerWidth;  // match container width dynamically
    const height = 400;

    // Create SVG with responsive sizing
    const svg = diagramDiv.append('svg')
        .attr('width', width)
        .attr('height', height)
        .attr('viewBox', `0 0 ${{width}} ${{height}}`)
        .style('overflow', 'visible')
        .style('display', 'block')
        .style('margin', '0 auto');
        
        // Get total reads for display
        const totalReads = Object.values(fusion.callers).filter(v => v !== null).reduce((a, b) => a + b, 0);
        
        // Encompassing reads indicator (top right)
        const readsGroup = svg.append('g').attr('transform', 'translate(1000, 60)');
        readsGroup.append('text')
            .attr('x', 0)
            .attr('y', 0)
            .style('font-size', '12px')
            .style('fill', '#2c3e50')
            .text(`Encompassing reads: ${{totalReads}}`);
        
        const exonWidth = 15;
        const exonHeight = 30;
        const exonGap = 2;
        
        /// === Draw Gene 1 (left side, right-aligned to breakpoint) ===
        const gene1Y = 200;
        let gene1X = 540;  // approximate breakpoint center alignment
        
        if (geneData.gene1 && geneData.gene1.exons) {{
            const g1 = geneData.gene1;
            const totalExons = g1.exons.length || 1;
            const bpIdx = g1.breakpoint_exon
                ? g1.exons.indexOf(g1.breakpoint_exon)
                : totalExons - 1;
        
            // Calculate total width up to breakpoint
            const exonWidth = 15;
            const exonGap = 2;
            const totalW = (bpIdx + 1) * (exonWidth + exonGap);
        
            // Shift gene1X left so that its breakpoint exon aligns under the fusion connector
            gene1X = 540 - totalW;
        }}
        
        if (geneData.gene1) {{
            drawGeneStructure(svg, geneData.gene1, gene1X, gene1Y, 15, 30, 2, '#3498db', 'left');
        }} else {{
            drawPlaceholderGene(svg, gene1, gene1X, gene1Y, 15, 30, 2, '#3498db', 'left');
        }}

        // Draw Gene 2 (right side) - Red
        const gene2X = 640;
        const gene2Y = 200;
        
        if (geneData.gene2) {{
            drawGeneStructure(svg, geneData.gene2, gene2X, gene2Y, exonWidth, exonHeight, exonGap, '#e74c3c', 'right');
        }} else {{
            drawPlaceholderGene(svg, gene2, gene2X, gene2Y, exonWidth, exonHeight, exonGap, '#e74c3c', 'right');
        }}
        // === Draw dashed line connecting breakpoint exons ===
        try {{
            if (geneData.gene1 && geneData.gene2) {{
                const g1 = geneData.gene1;
                const g2 = geneData.gene2;
        
                const exons1 = g1.exons || [];
                const exons2 = g2.exons || [];
        
                // Get indices of breakpoint exons
                const bpIdx1 = g1.breakpoint_exon ? exons1.indexOf(g1.breakpoint_exon) : exons1.length - 1;
                const bpIdx2 = g2.breakpoint_exon ? exons2.indexOf(g2.breakpoint_exon) : 0;
        
                // Constants for exon box layout
                const exonW = 15, exonGap = 2;
        
                // --- Calculate X positions at exon *edges* ---
                const bpX1 = gene1X + (bpIdx1 + 1) * (exonW + exonGap);  // right edge of gene1 breakpoint exon
                const bpX2 = gene2X;                                     // left edge of gene2 breakpoint exon
        
                // Y coordinates for connection (at bottom of exon boxes)
                const bpY = gene1Y + exonHeight / 2;
        
                // Draw dashed connector line
                svg.append('line')
                    .attr('x1', bpX1)
                    .attr('y1', bpY)
                    .attr('x2', bpX2)
                    .attr('y2', bpY)
                    .attr('stroke', '#2c3e50')
                    .attr('stroke-width', 2)
                    .attr('stroke-dasharray', '6,4');
            }}
        }} catch (err) {{
            console.error('Error drawing dashed connector:', err);
        }}

        // Draw fusion connection at top (simplified; avoids undefined variables)
        drawFusionConnection(svg, geneData, exonWidth, exonHeight);
        
        // === Draw mini cytoband chromosomes with breakpoints ===
        drawCytobandMiniChromosomes(svg, geneData, width, height);
}}

// === Transcript dropdowns for both genes ===
let currentFusion = null;
let selectedTranscripts = {{ gene1: null, gene2: null }};

// === Populate transcript dropdowns dynamically ===
function populateTranscriptDropdowns(fusion) {{
    const gene1Select = document.getElementById('gene1TranscriptSelect');
    const gene2Select = document.getElementById('gene2TranscriptSelect');

    const fusionName = fusion.fusion;
    const gStruct = geneStructureData[fusionName];
    if (!gStruct) {{
        console.warn(`No gene structure found for fusion: ${{fusionName}}`);
        return;
    }}

    // Helper to fill a dropdown
    function fillSelect(selectElem, geneData, color) {{
        selectElem.innerHTML = '';
        if (!geneData || !geneData.available_transcripts || !geneData.available_transcripts.length) {{
            selectElem.disabled = true;
            const opt = document.createElement('option');
            opt.textContent = 'N/A';
            selectElem.appendChild(opt);
            return;
        }}

        selectElem.disabled = false;

        geneData.available_transcripts.forEach(tid => {{
            const opt = document.createElement('option');
            opt.value = tid;
            opt.textContent = tid;
            if (tid === geneData.transcript) opt.selected = true;
            selectElem.appendChild(opt);
        }});

        // Highlight default transcript visually
        selectElem.style.border = `2px solid ${{color}}`;
        selectElem.style.fontWeight = '600';
    }}

    fillSelect(gene1Select, gStruct.gene1, '#667eea'); // blue for gene1
    fillSelect(gene2Select, gStruct.gene2, '#e74c3c'); // red for gene2
}}

// === When user changes transcript ===
function updateFusionTranscript(side, newTranscriptId) {{
    if (!selectedFusion) return;
    const fusionName = selectedFusion;
    const gStruct = geneStructureData[fusionName];
    if (!gStruct || !gStruct[side]) return;

    const geneData = gStruct[side];
    const oldTranscript = geneData.transcript;
    geneData.transcript = newTranscriptId;

    // --- Load exon data for this transcript ---
    if (geneData.all_transcripts && geneData.all_transcripts[newTranscriptId]) {{
        const tdata = geneData.all_transcripts[newTranscriptId];
        const exonData = tdata.exons || [];
        geneData.exons = exonData.map(e => e.exon_number);
        geneData.total_exons = geneData.exons.length;
        
        // --- CRITICAL: Recalculate breakpoint exon for the new transcript ---
        const breakpointPos = geneData.pos;
        let bpExon = null;
        
        // Find which exon contains the breakpoint position
        for (const exon of exonData) {{
            if (exon.start <= breakpointPos && breakpointPos <= exon.end) {{
                bpExon = exon.exon_number;
                break;
            }}
        }}
        
        // If not inside any exon, find the closest exon
        if (bpExon === null && exonData.length > 0) {{
            const closest = exonData.reduce((prev, curr) => {{
                const prevDist = Math.min(Math.abs(prev.start - breakpointPos), Math.abs(prev.end - breakpointPos));
                const currDist = Math.min(Math.abs(curr.start - breakpointPos), Math.abs(curr.end - breakpointPos));
                return currDist < prevDist ? curr : prev;
            }});
            bpExon = closest.exon_number;
        }}
        
        geneData.breakpoint_exon = bpExon;
        console.log(`Recalculated breakpoint exon for ${{side}}: exon ${{bpExon}} at position ${{breakpointPos}}`);
    }} else {{
        console.warn(`No exon data found for transcript ${{newTranscriptId}} of ${{side}}`);
    }}

    // --- Redraw the gene fusion diagram with the new exon structure ---
    const fusion = fusionData.find(f => f.fusion === fusionName);
    if (fusion) {{
        const [gene1, gene2] = fusion.fusion.split(':');
        drawGeneFusionDiagram(fusion, gene1, gene2);
    }}

    console.log(`Updated ${{side}} transcript: ${{oldTranscript}} ${{newTranscriptId}} (breakpoint exon: ${{geneData.breakpoint_exon}})`);
}}


// === Simple fallback for missing gene data ===
function drawPlaceholderGene(svg, x, y, color = '#bdc3c7') {{
    const width = 80;
    const height = 20;
    svg.append('rect')
        .attr('x', x)
        .attr('y', y)
        .attr('width', width)
        .attr('height', height)
        .attr('fill', color)
        .attr('opacity', 0.3)
        .attr('rx', 4)
        .attr('ry', 4);
    svg.append('text')
        .attr('x', x + width / 2)
        .attr('y', y + height / 1.5)
        .attr('text-anchor', 'middle')
        .style('font-size', '11px')
        .style('fill', '#2c3e50')
        .text('No Data');
}}


function getCytobandForPosition(chr, pos) {{
    if (!Array.isArray(cytobandData) || cytobandData.length === 0) return 'N/A';
    
    // Cytoband file uses chr-prefixed names (e.g. chr1)
    const chrName = chr.startsWith('chr') ? chr : `chr${{chr}}`;
    
    // Find the matching cytoband
    const band = cytobandData.find(b => 
        b.chrom === chrName && pos >= b.start && pos <= b.end
    );
    
    return band ? band.name : 'N/A';
}}

function drawFusionConnection(svg, geneData, exonWidth, exonHeight) {{
    if (!geneData || !geneData.gene1 || !geneData.gene2) return;

    const gene1 = geneData.gene1;
    const gene2 = geneData.gene2;

    const exons1 = gene1.exons || [];
    const exons2 = gene2.exons || [];

    // === Handle missing breakpoint exons robustly ===
    let bpIdx1 = exons1.indexOf(gene1.breakpoint_exon);
    if (bpIdx1 === -1 || bpIdx1 === undefined) {{
        // fallback to middle exon if missing
        bpIdx1 = Math.floor(exons1.length / 2);
        console.warn(`Fallback: breakpoint exon not found in gene1 (${{gene1.name || 'unknown'}}), using exon index ${{bpIdx1}}`);
    }}

    let bpIdx2 = exons2.indexOf(gene2.breakpoint_exon);
    if (bpIdx2 === -1 || bpIdx2 === undefined) {{
        bpIdx2 = Math.floor(exons2.length / 2);
        console.warn(`Fallback: breakpoint exon not found in gene2 (${{gene2.name || 'unknown'}}), using exon index ${{bpIdx2}}`);
    }}

    // === Simplified schematic (top, context-aware + centered) ===
    const exonW = 20, exonGap = 4, exonH = 30;
    const schematicY = 50;

    // ----- Left gene -----
    let leftBlocks = [];
    if (bpIdx1 > 0) {{
        leftBlocks = [
            {{ label: "1", color: "#3498db" }},
            {{ label: "...", color: "none" }},
            {{ label: gene1.breakpoint_exon ? gene1.breakpoint_exon.toString() : "BP", color: "#3498db" }}
        ];
    }} else {{
        leftBlocks = [
            {{ label: gene1.breakpoint_exon ? gene1.breakpoint_exon.toString() : "BP", color: "#3498db" }}
        ];
    }}

    // ----- Right gene -----
    let rightBlocks = [];
    if (bpIdx2 < exons2.length - 1) {{
        rightBlocks = [
            {{ label: gene2.breakpoint_exon ? gene2.breakpoint_exon.toString() : "BP", color: "#e74c3c" }},
            {{ label: "...", color: "none" }},
            {{ label: exons2.length.toString(), color: "#e74c3c" }}
        ];
    }} else {{
        rightBlocks = [
            {{ label: gene2.breakpoint_exon ? gene2.breakpoint_exon.toString() : "BP", color: "#e74c3c" }}
        ];
    }}

    // --- Center schematic horizontally ---
    const totalBlocks = leftBlocks.length + rightBlocks.length + 1;
    const schematicWidth = totalBlocks * (exonW + exonGap);
    const svgWidth = parseFloat(svg.attr('width')) || 1200;
    const startX = (svgWidth - schematicWidth) / 2;

    // ----- Draw left schematic -----
    leftBlocks.forEach((block, i) => {{
        const x = startX + i * (exonW + exonGap);
        if (block.color === "none") {{
            svg.append('text')
                .attr('x', x + exonW / 2)
                .attr('y', schematicY + 18)
                .attr('text-anchor', 'middle')
                .style('font-weight', 'bold')
                .text('...');
        }} else {{
            svg.append('rect')
                .attr('x', x).attr('y', schematicY)
                .attr('width', exonW).attr('height', exonH)
                .attr('fill', block.color).attr('stroke', 'white').attr('stroke-width', 2);
            svg.append('text')
                .attr('x', x + exonW / 2).attr('y', schematicY + 17)
                .attr('text-anchor', 'middle')
                .style('fill', 'white').style('font-weight', 'bold')
                .text(block.label);
        }}
    }});

    // ----- Draw right schematic -----
    const rightStart = startX + (leftBlocks.length + 1) * (exonW + exonGap);
    rightBlocks.forEach((block, i) => {{
        const x = rightStart + i * (exonW + exonGap);
        if (block.color === "none") {{
            svg.append('text')
                .attr('x', x + exonW / 2)
                .attr('y', schematicY + 18)
                .attr('text-anchor', 'middle')
                .style('font-weight', 'bold')
                .text('...');
        }} else {{
            svg.append('rect')
                .attr('x', x).attr('y', schematicY)
                .attr('width', exonW).attr('height', exonH)
                .attr('fill', block.color).attr('stroke', 'white').attr('stroke-width', 2);
            svg.append('text')
                .attr('x', x + exonW / 2).attr('y', schematicY + 17)
                .attr('text-anchor', 'middle')
                .style('fill', 'white').style('font-weight', 'bold')
                .text(block.label);
        }}
    }});

    // Label
    svg.append('text')
        .attr('x', startX - 60)
        .attr('y', schematicY + 18)
        .attr('text-anchor', 'end')
        .style('font-size', '14px')
        .style('font-weight', 'bold')
        .style('fill', '#2c3e50')
        .text('Fusion:');

    // Dashed connector down to detailed structures
    const bpLeftX = startX + Math.max(0, leftBlocks.length - 1) * (exonW + exonGap) + exonW / 2;
    const bpRightX = rightStart + exonW / 2;
    const fuseMidX = (bpLeftX + bpRightX) / 2;

    if (isNaN(fuseMidX)) {{
        console.warn("Skipping fusion connector: invalid coordinates", bpLeftX, bpRightX);
        return;
    }}

    svg.append('line')
        .attr('x1', fuseMidX)
        .attr('y1', schematicY + exonH)
        .attr('x2', fuseMidX)
        .attr('y2', schematicY + 80)
        .attr('stroke', '#2c3e50')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '5,5');
}}

function drawGeneStructure(svg, geneData, x, y, exonWidth, exonHeight, exonGap, color, side) {{
    const group = svg.append('g').attr('transform', `translate(${{x}}, ${{y}})`);

    const allExons = geneData.exons || [];
    const bpIdx = geneData.breakpoint_exon !== undefined && geneData.breakpoint_exon !== null
        ? allExons.indexOf(geneData.breakpoint_exon)
        : (side === 'left' ? allExons.length - 1 : 0);

    // Select subset for fusion region
    let exons = [];
    if (side === 'left') {{
        exons = allExons.slice(0, bpIdx + 1);
    }} else {{
        exons = allExons.slice(bpIdx);
    }}

    // Gene line
    const totalW = exons.length * (exonWidth + exonGap);
    group.append('line')
        .attr('x1', 0).attr('y1', exonHeight / 2)
        .attr('x2', totalW)
        .attr('y2', exonHeight / 2)
        .attr('stroke', '#2c3e50').attr('stroke-width', 2);

    // Draw retained exons
    exons.forEach((exon, i) => {{
        const isBP = (side === 'left' && i === exons.length - 1) ||
                     (side === 'right' && i === 0);

        group.append('rect')
            .attr('x', i * (exonWidth + exonGap))
            .attr('y', 0)
            .attr('width', exonWidth)
            .attr('height', exonHeight)
            .attr('fill', color)
            .attr('opacity', isBP ? 1 : 0.95)
            .attr('stroke', isBP ? '#2c3e50' : 'white')
            .attr('stroke-width', isBP ? 3 : 1.5);

        group.append('text')
            .attr('x', i * (exonWidth + exonGap) + exonWidth / 2)
            .attr('y', exonHeight / 2 + 4)
            .attr('text-anchor', 'middle')
            .style('fill', 'white')
            .style('font-weight', 'bold')
            .style('font-size', '10px')
            .text(exon);
    }});

    // Dashed breakpoint connector
    const bpX = (side === 'left')
        ? (exons.length - 1) * (exonWidth + exonGap) + exonWidth / 2
        : exonWidth / 2;

    group.append('line')
        .attr('x1', bpX).attr('y1', exonHeight)
        .attr('x2', bpX).attr('y2', exonHeight + 40)
        .attr('stroke', color).attr('stroke-width', 2)
        .attr('stroke-dasharray', '5,5');

    // === Breakpoint coordinate (above gene) ===
    const bpTextY = y - 15;
    const chr = geneData.chr || 'chr?';
    const pos = geneData.pos || 0;

    svg.append('text')
        .attr('x', x + totalW / 2)
        .attr('y', bpTextY)
        .attr('text-anchor', 'middle')
        .style('font-size', '11px')
        .style('fill', color)
        .style('font-weight', 'bold')
        .text(`${{chr.replace('chr', '')}}:${{pos.toLocaleString()}}`);

    // === Labels below ===
    svg.append('text')
        .attr('x', x + totalW / 2)
        .attr('y', y + exonHeight + 70)
        .attr('text-anchor', 'middle')
        .style('font-size', '14px')
        .style('font-weight', 'bold')
        .style('fill', '#2c3e50')
        .text(geneData.name);

    svg.append('text')
        .attr('x', x + totalW / 2)
        .attr('y', y + exonHeight + 88)
        .attr('text-anchor', 'middle')
        .style('font-size', '11px')
        .style('fill', '#2c3e50')
        .style('font-weight', '600')
        .text(geneData.transcript || '');

    // === Cytoband (bottom, below transcript) ===
    //const cytoband = getCytobandForPosition(chr, pos);

    //svg.append('text')
        //.attr('x', x + totalW / 2)
        //.attr('y', y + exonHeight + 105)
        //.attr('text-anchor', 'middle')
        //.style('font-size', '11px')
        //.style('fill', '#000000')
        //.style('font-weight', 'bold')
        //.text(cytoband && cytoband !== 'N/A'
            //? `${{chr}}${{cytoband.replace(/^chr/, '')}}`
            //: '');
}}

function drawCytobandMiniChromosomes(svg, geneData, width, height) {{
    const chromHeight = 12;
    const chromWidth = 120;
    const verticalOffset = 360;  // position above exon structures

    const genes = [
        {{ gene: geneData.gene1, color: '#3498db', side: 'left' }},
        {{ gene: geneData.gene2, color: '#e74c3c', side: 'right' }}
    ].filter(d => d.gene && d.gene.chr);

    genes.forEach(({{ gene, color, side }}) => {{
        const chrName = gene.chr.startsWith('chr') ? gene.chr : `chr${{gene.chr}}`;
        const chrBands = cytobandData.filter(b => b.chrom === chrName);
        if (!chrBands.length) return;

        const chrStart = Math.min(...chrBands.map(b => b.start));
        const chrEnd = Math.max(...chrBands.map(b => b.end));
        const scale = d3.scaleLinear().domain([chrStart, chrEnd]).range([0, chromWidth]);

        const geneCenterX =
            side === 'left' ? 540 - chromWidth / 2 - 30 : 640 - chromWidth / 2 + 30;
        const yBase = verticalOffset;

        // === Draw cytoband rectangles ===
        chrBands.forEach(band => {{
            const x = geneCenterX + scale(band.start);
            const w = scale(band.end) - scale(band.start);
            const fill = band.stain.includes('neg') ? '#ffffff'
                : band.stain.includes('pos25') ? '#d0d0d0'
                : band.stain.includes('pos50') ? '#a0a0a0'
                : band.stain.includes('pos75') ? '#707070'
                : band.stain.includes('pos100') ? '#404040'
                : band.stain.includes('acen') ? '#d66'
                : '#e0e0e0';

            svg.append('rect')
                .attr('x', x)
                .attr('y', yBase)
                .attr('width', w)
                .attr('height', chromHeight)
                .attr('fill', fill)
                .attr('stroke', '#555')   // darker boundary
                .attr('stroke-width', 1.2); // thicker outline
        }});

        // === Telomere semicircles (ends of chromosome) ===
        const radius = chromHeight / 2;
        const leftCenterX = geneCenterX;
        const rightCenterX = geneCenterX + chromWidth;
        const centerY = yBase + radius;

        // Left telomere
        svg.append('path')
            .attr('d', `M ${{leftCenterX}} ${{yBase}}
                        A ${{radius}} ${{radius}} 0 0 0 ${{leftCenterX}} ${{yBase + chromHeight}} Z`)
            .attr('fill', '#ccc')
            .attr('stroke', '#555')
            .attr('stroke-width', 1.2);

        // Right telomere
        svg.append('path')
            .attr('d', `M ${{rightCenterX}} ${{yBase}}
                        A ${{radius}} ${{radius}} 0 0 1 ${{rightCenterX}} ${{yBase + chromHeight}} Z`)
            .attr('fill', '#ccc')
            .attr('stroke', '#555')
            .attr('stroke-width', 1.2);

        // === Centromere triangles ===
        const acenBands = chrBands.filter(b => b.stain && b.stain.includes('acen'));
        if (acenBands.length >= 2) {{
            const cenStart = geneCenterX + scale(Math.min(...acenBands.map(b => b.start)));
            const cenEnd = geneCenterX + scale(Math.max(...acenBands.map(b => b.end)));

            const mid = (cenStart + cenEnd) / 2;
            const topY = yBase;
            const bottomY = yBase + chromHeight;

            // Left triangle
            svg.append('path')
                .attr('d', `M ${{cenStart}} ${{topY}} L ${{mid}} ${{centerY}} L ${{cenStart}} ${{bottomY}} Z`)
                .attr('fill', '#b44')
                .attr('stroke', '#555')
                .attr('stroke-width', 0.8);

            // Right triangle
            svg.append('path')
                .attr('d', `M ${{cenEnd}} ${{topY}} L ${{mid}} ${{centerY}} L ${{cenEnd}} ${{bottomY}} Z`)
                .attr('fill', '#b44')
                .attr('stroke', '#555')
                .attr('stroke-width', 0.8);
        }}

        // === Breakpoint tick (solid red) ===
        const pos = gene.pos || (chrStart + chrEnd) / 2;
        const bpX = geneCenterX + scale(pos);
        svg.append('line')
            .attr('x1', bpX)
            .attr('x2', bpX)
            .attr('y1', yBase - 5)
            .attr('y2', yBase + chromHeight + 5)
            .attr('stroke', '#e74c3c')
            .attr('stroke-width', 2);

        // === Cytoband label (below bar) ===
        const bandName = getCytobandForPosition(chrName, pos);
        svg.append('text')
            .attr('x', geneCenterX + chromWidth / 2)
            .attr('y', yBase + chromHeight + 18)
            .attr('text-anchor', 'middle')
            .style('font-size', '11px')
            .style('font-weight', 'bold')
            .style('fill', '#2c3e50')
            .text(`${{chrName.replace('chr', '')}}${{bandName ? bandName.replace(/^chr/, '') : ''}}`);
    }});
}}

function drawCoveragePlot(data, gene1, gene2) {{
    const plotDiv = d3.select('#coveragePlot');
    plotDiv.selectAll('*').remove();
    
    const margin = {{top: 20, right: 30, bottom: 100, left: 60}};
    const width = 800 - margin.left - margin.right;
    const height = 400 - margin.top - margin.bottom;
    
    const svg = plotDiv.append('svg')
        .attr('width', width + margin.left + margin.right)
        .attr('height', height + margin.top + margin.bottom)
        .append('g')
        .attr('transform', `translate(${{margin.left}},${{margin.top}})`);
    
    const x = d3.scaleBand()
        .domain(data.map(d => d.region))
        .range([0, width])
        .padding(0.2);
    
    const y = d3.scaleLinear()
        .domain([0, d3.max(data, d => parseFloat(d.mean_coverage || d.mean || 0))])
        .nice()
        .range([height, 0]);
    
    svg.selectAll('.bar')
        .data(data)
        .enter()
        .append('rect')
        .attr('class', 'bar')
        .attr('x', d => x(d.region))
        .attr('width', x.bandwidth())
        .attr('y', d => y(parseFloat(d.mean_coverage || d.mean || 0)))
        .attr('height', d => height - y(parseFloat(d.mean_coverage || d.mean || 0)))
        .attr('fill', d => d.region.includes(gene1) ? '#667eea' : '#764ba2');
    
    svg.append('g')
        .attr('transform', `translate(0,${{height}})`)
        .call(d3.axisBottom(x))
        .selectAll('text')
        .attr('transform', 'rotate(-45)')
        .style('text-anchor', 'end');
    
    svg.append('g')
        .call(d3.axisLeft(y));
    
    svg.append('text')
        .attr('x', width / 2)
        .attr('y', -5)
        .attr('text-anchor', 'middle')
        .style('font-size', '14px')
        .style('font-weight', 'bold')
        .text('Mean Coverage by Exon');
    
    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('y', -40)
        .attr('x', -height / 2)
        .attr('text-anchor', 'middle')
        .text('Mean Coverage');
}}

function switchCoverageTab(tab, evt) {{
    // Update tab buttons
    document.querySelectorAll('.coverage-tab').forEach(btn => {{
        btn.classList.remove('active');
    }});
    if(evt && evt.target) evt.target.classList.add('active');
    
    // Show/hide tables
    if (tab === 'mosdepth') {{
        document.getElementById('mosdepthTable').style.display = 'block';
        document.getElementById('bedtoolsTable').style.display = 'none';
    }} else {{
        document.getElementById('mosdepthTable').style.display = 'none';
        document.getElementById('bedtoolsTable').style.display = 'block';
    }}
}}

function drawMosdepthTable(data) {{
    const tableDiv = d3.select('#mosdepthTable');
    tableDiv.selectAll('*').remove();
    
    if (!data || data.length === 0) {{
        tableDiv.append('p')
            .style('text-align', 'center')
            .style('padding', '20px')
            .style('color', '#7f8c8d')
            .text('No mosdepth coverage data available');
        return;
    }}
    
    // Create table structure
    const table = tableDiv.append('table')
        .style('width', '100%')
        .style('border-collapse', 'collapse')
        .style('margin-top', '10px');
    
    // Get column names - columns 1-4 and 11 (0-indexed: 0-3 and 10)
    const firstRow = data[0];
    const allKeys = Object.keys(firstRow);
    
    // Select columns: first 4 and column 11 (index 10)
    const selectedKeys = allKeys.slice(0, 4);
    if (allKeys.length > 10) {{
        selectedKeys.push(allKeys[10]);
    }}
    
    // Header
    const thead = table.append('thead');
    const headerRow = thead.append('tr');
    selectedKeys.forEach(key => {{
        headerRow.append('th')
            .style('background-color', '#667eea')
            .style('color', 'white')
            .style('padding', '10px')
            .style('text-align', 'left')
            .style('font-weight', '600')
            .text(key);
    }});
    
    // Body
    const tbody = table.append('tbody');
    data.forEach((row, i) => {{
        const tr = tbody.append('tr');
        if (i % 2 === 0) {{
            tr.style('background-color', '#f8f9fa');
        }}
        
        selectedKeys.forEach(key => {{
            tr.append('td')
                .style('padding', '8px 10px')
                .style('border-bottom', '1px solid #e0e0e0')
                .text(row[key] || '-');
        }});
    }});
}}

function drawBedtoolsTable(data) {{
    const tableDiv = d3.select('#bedtoolsTable');
    tableDiv.selectAll('*').remove();
    
    if (!data || data.length === 0) {{
        tableDiv.append('p')
            .style('text-align', 'center')
            .style('padding', '20px')
            .style('color', '#7f8c8d')
            .text('No bedtools coverage data available');
        return;
    }}
    
    // Create table structure
    const table = tableDiv.append('table')
        .style('width', '100%')
        .style('border-collapse', 'collapse')
        .style('margin-top', '10px');
    
    // Get ALL column names
    const firstRow = data[0];
    const allKeys = Object.keys(firstRow);
    
    // Header
    const thead = table.append('thead');
    const headerRow = thead.append('tr');
    allKeys.forEach(key => {{
        headerRow.append('th')
            .style('background-color', '#764ba2')
            .style('color', 'white')
            .style('padding', '10px')
            .style('text-align', 'left')
            .style('font-weight', '600')
            .text(key);
    }});
    
    // Body
    const tbody = table.append('tbody');
    data.forEach((row, i) => {{
        const tr = tbody.append('tr');
        if (i % 2 === 0) {{
            tr.style('background-color', '#f8f9fa');
        }}
        
        allKeys.forEach(key => {{
            tr.append('td')
                .style('padding', '8px 10px')
                .style('border-bottom', '1px solid #e0e0e0')
                .text(row[key] || '-');
        }});
    }});
}}

function closeDetails() {{
    selectedFusion = null;
    
    // Reset circos highlighting
    d3.selectAll('.fusion-arc').attr('opacity', 0.6).attr('stroke-width', function() {{
        const currentWidth = parseFloat(d3.select(this).attr('stroke-width'));
        return Math.max(1.5, currentWidth - 2);
    }});
    document.getElementById('detailsPanel').classList.remove('active');
    document.getElementById('coverageSection').style.display = 'none';
}}

// Wire up controls
document.getElementById('minCallers').addEventListener('change', updateTable);
document.getElementById('minReads').addEventListener('change', updateTable);
document.getElementById('sortBy').addEventListener('change', updateTable);

// Initial render
updateTable();
</script>
</body>
</html>
"""

with open(output_file,  'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"Wrote HTML dashboard to {output_file}")
