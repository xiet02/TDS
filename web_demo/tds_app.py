"""
TDS: Target Discovery System - Oncology Target Discovery Chatbot
==================================================================

Companion script for the TAA paper:
"Target Discovery System: A GraphRAG-based approach for oncology target discovery"
https://www.biorxiv.org/content/10.1101/2025.05.06.652559v1

This Gradio-based web interface provides an interactive chatbot for querying
the GraphRAG knowledge graph built from oncology literature. Users can ask
natural language questions about cancer targets, safety profiles, and 
therapeutic potential.

Features:
---------
- Interactive chat interface powered by GraphRAG
- Support for both local and global search methods
- Citation tracking and source visualization
- Real-time query processing with streaming output
- Dark-themed UI optimized for biomedical research

Requirements:
-------------
- gradio
- pandas
- graphrag (Microsoft GraphRAG)
- conda environment named 'graphrag_env'

Usage:
------
    python tds_chatbot.py

Environment Variables:
----------------------
- GRAPH_RAG_ROOT: Path to GraphRAG root directory
- GRAPHRAG_PYTHON: Path to Python executable in graphrag_env
- CONDA_EXE: Path to conda executable
- SHOW_THINKING: Set to "1" to display query execution details
- SHOW_SOURCES: Set to "1" to display artifact citations

Author: xiet02
License: MIT
Repository: https://github.com/xiet02/TDS
"""

import gradio as gr
import subprocess
import os
import re
import pandas as pd
from pathlib import Path

# ============================================================================
# Citation and Data Loading Functions
# ============================================================================

def _parse_citations(text):
    """
    Extract parquet file citations from GraphRAG output.
    
    GraphRAG generates intermediate artifacts (parquet files) during query
    processing. This function identifies references to these artifacts in the
    output text for source tracking.
    
    Args:
        text (str): Raw output text from GraphRAG query
        
    Returns:
        list: Unique list of parquet file paths found in the text
    """
    # Pattern 1: Unix-style paths to output artifacts
    # Pattern 2: Windows-style paths to output artifacts
    patterns = [
        r"output/artifacts/[\w\-/]+\.parquet",
        r"[A-Z]:\\[^\n]*?output\\artifacts\\[^\n]*?\.parquet",
    ]
    
    found = []
    for pattern in patterns:
        found.extend(re.findall(pattern, text))
    
    # Remove duplicates while preserving order
    return list(dict.fromkeys(found))


def _load_parquet(path):
    """
    Safely load a parquet file with error handling.
    
    Args:
        path (str): File path to parquet file
        
    Returns:
        pd.DataFrame or None: DataFrame if successful, None otherwise
    """
    if not os.path.exists(path):
        return None
    
    try:
        return pd.read_parquet(path)
    except Exception:
        return None


# ============================================================================
# Configuration and Environment Setup
# ============================================================================

# GraphRAG root directory - contains the indexed knowledge graph
GRAPH_RAG_ROOT = os.environ.get("GRAPH_RAG_ROOT", r"path-to-graphrag")

# Logo for the TDS interface
LOGO_PATH = str(Path(__file__).parent / "assets" / "tds_logo.svg")

# Python executable within the graphrag conda environment
GRAPHRAG_PYTHON = os.environ.get("GRAPHRAG_PYTHON", r"path-to-python.exe")

# Fallback paths for conda executable (platform-specific)
CONDA_EXE_DEFAULTS = [
    r"path-to-python.exe",  # Update with your actual default path
]

# Attempt to locate conda executable
CONDA_EXE = os.environ.get("CONDA_EXE", "")
if not CONDA_EXE or not os.path.exists(CONDA_EXE):
    for _p in CONDA_EXE_DEFAULTS:
        if os.path.exists(_p):
            CONDA_EXE = _p
            break

# Feature flags for debugging and transparency
SHOW_THINKING = os.environ.get("SHOW_THINKING", "0") == "1"  # Show query execution details
SHOW_SOURCES = os.environ.get("SHOW_SOURCES", "0") == "1"    # Show artifact citations


# ============================================================================
# Utility Functions
# ============================================================================

def _fmt_cmd(cmd_list):
    """
    Format command list for display, quoting arguments with spaces.
    
    Args:
        cmd_list (list): List of command components
        
    Returns:
        str: Formatted command string
    """
    def q(x):
        if any(ch.isspace() for ch in x):
            return f'"{x}"'
        return x
    
    return " ".join(q(str(s)) for s in cmd_list)


def _build_df(message, output_text):
    """
    Build a DataFrame of relevant nodes and entities from GraphRAG artifacts.
    
    This function searches the GraphRAG knowledge graph for nodes and entities
    that match tokens found in the user's query and the response. It's useful
    for displaying relevant context and provenance.
    
    Args:
        message (str): User's query text
        output_text (str): GraphRAG response text
        
    Returns:
        pd.DataFrame: Combined DataFrame of matching nodes and entities
    """
    # Define paths to GraphRAG artifact files
    base = os.path.join(GRAPH_RAG_ROOT, "output", "artifacts")
    nodes_path = os.path.join(base, "create_final_nodes.parquet")
    entities_path = os.path.join(base, "create_final_entities.parquet")
    
    # Extract tokens (3+ alphanumeric characters) from query and response
    tokens = set(re.findall(r"[A-Za-z0-9]{3,}", f"{message} {output_text}"))
    
    # Build regex pattern for token matching (longest tokens first for better matching)
    if tokens:
        regex = "|".join(sorted({re.escape(t) for t in tokens}, key=len, reverse=True))
    else:
        regex = None
    
    frames = []
    
    # Process both nodes and entities
    for label, p in [("nodes", nodes_path), ("entities", entities_path)]:
        df = _load_parquet(p)
        if df is None:
            continue
        
        # Filter DataFrame based on token matches in string columns
        if regex:
            mask = False
            for c in df.columns:
                if df[c].dtype == object:
                    m = df[c].astype(str).str.contains(regex, case=False, na=False)
                    mask = m if isinstance(mask, bool) else (mask | m)
            
            if isinstance(mask, bool):
                filtered = df.head(10)  # No matches, show first 10
            else:
                filtered = df[mask].head(20)  # Show matching rows
        else:
            filtered = df.head(10)  # No tokens, show first 10
        
        # Add source label for tracking
        filtered = filtered.copy()
        filtered["source"] = label
        frames.append(filtered)
    
    # Combine all DataFrames
    if not frames:
        return pd.DataFrame()
    
    return pd.concat(frames, ignore_index=True)


# ============================================================================
# Main Query Function
# ============================================================================

def graphrag_query(message, history, method):
    """
    Execute a GraphRAG query and return formatted results.
    
    This is the core function that interfaces with the GraphRAG CLI tool.
    It handles query execution, output parsing, citation extraction, and
    error handling.
    
    Args:
        message (str): User's query text
        history (list): Conversation history (unused in current implementation)
        method (str): Search method ('local' or 'global')
        
    Returns:
        tuple: (response_text, overlay_update, modal_update, button_update)
            - response_text: Formatted response with citations and details
            - overlay_update: Gradio update for error overlay
            - modal_update: Gradio update for error modal
            - button_update: Gradio update for close button
    """
    # Currently hardcoded to global search (can be parameterized if needed)
    method_used = "global"
    acc = ""  # Accumulator for stdout
    
    # Build GraphRAG query command
    cmd_main = [
        CONDA_EXE if CONDA_EXE else "conda",  # Conda executable
        "run", "-n", "graphrag_env",           # Run in graphrag_env
        "--no-capture-output",                 # Stream output
        "graphrag", "query",                   # GraphRAG query command
        "--root", GRAPH_RAG_ROOT,              # Knowledge graph root
        "--method", method_used,               # Search method
        "--query", message,                    # User query
    ]
    
    # Execute GraphRAG query
    try:
        proc = subprocess.Popen(
            cmd_main,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False
        )
        cmd_used = _fmt_cmd(cmd_main)
    except Exception as e:
        # Handle execution errors
        err_text = str(e)
        modal = f"### Query Error\n\n{err_text}\n\nCommand: {_fmt_cmd(cmd_main)}"
        return (
            None,
            gr.update(value="", visible=True),    # Show overlay
            gr.update(value=modal, visible=True),  # Show error modal
            gr.update(visible=True)                # Show close button
        )
    
    # Stream stdout line by line
    if proc.stdout is not None:
        for line in proc.stdout:
            acc += line
    
    # Wait for process completion
    rc = proc.wait()
    err = proc.stderr.read() if proc.stderr is not None else ""
    
    # Handle non-zero return codes
    if rc != 0:
        if not acc.strip():
            modal = f"### Query Error\n\n{(err or 'GraphRAG query failed').strip()}\n\nCommand: {cmd_used}"
            return (
                None,
                gr.update(value="", visible=True),
                gr.update(value=modal, visible=True),
                gr.update(visible=True)
            )
    
    # Parse citations from output
    citations = _parse_citations(acc)
    
    # Build thinking/debug section
    thinking = f"Command: {cmd_used}"
    
    # Build sources section
    sources_md = "\n".join([f"- {c}" for c in citations]) if citations else "- No artifact references found."
    
    # Assemble optional details sections
    details = ""
    
    if SHOW_THINKING:
        details += "\n\n<details><summary>Thinking</summary>\n" + thinking + "\n</details>\n"
    
    if SHOW_SOURCES:
        details += "<details><summary>Sources</summary>\n" + sources_md + "\n</details>\n"
    
    # Include warnings/errors if present
    if err.strip():
        details += (
            "\n\n<details><summary>Warnings</summary>\n" + err.strip() + "\n</details>\n"
        )
    
    # Format final response
    reply_text = acc.strip() if acc.strip() else ""
    final_msg = (reply_text + details) if reply_text else details
    
    return (
        final_msg,
        gr.update(visible=False),  # Hide overlay
        gr.update(visible=False),  # Hide modal
        gr.update(visible=False)   # Hide close button
    )


# ============================================================================
# Gradio Interface
# ============================================================================

# Custom theme for the interface
custom_theme = gr.themes.Soft()

with gr.Blocks() as demo:
    # ========================================================================
    # Custom CSS and JavaScript
    # ========================================================================
    gr.HTML("""
    <style>
    /* Hide unnecessary Gradio controls */
    button[aria-label="Download"],
    button[aria-label="Fullscreen"],
    button[aria-label="Zoom in"],
    button[aria-label="Zoom out"],
    button[aria-label="Share"] { display: none !important; }
    
    /* Error modal styling */
    #tds-error-overlay { 
        position: fixed; 
        inset: 0; 
        background: rgba(0,0,0,0.35); 
        z-index: 999; 
    }
    #tds-error-modal { 
        position: fixed; 
        top: 50%; 
        left: 50%; 
        transform: translate(-50%, -50%); 
        background: #ffffff; 
        border: 2px solid #0f9d8f; 
        border-radius: 12px; 
        padding: 16px; 
        z-index: 1000; 
        box-shadow: 0 10px 30px rgba(0,0,0,0.25); 
        max-width: 640px; 
    }
    #tds-error-modal h3 { margin-top: 0; color: #d32f2f; }
    #tds-error-modal h4 { margin-top: 0; color: #0f9d8f; }
    #tds-error-close { 
        position: fixed; 
        top: calc(50% - 64px); 
        left: calc(50% + 292px); 
        transform: translate(-50%, -50%); 
        z-index: 1001; 
        background: #0f9d8f; 
        color: #fff; 
        border-radius: 50%; 
        width: 28px; 
        height: 28px; 
        line-height: 24px; 
        text-align: center; 
        padding: 0; 
    }
    
    /* Header styling */
    #tds-header { 
        display: flex; 
        align-items: center; 
        gap: 8px; 
    }
    
    /* Dark theme for biomedical research aesthetic */
    body, .gradio-container { 
        background: #0e1217 !important; 
        color: #e6edf3 !important; 
    }
    .gradio-container .block, 
    .gradio-container .panel, 
    .gradio-container .wrap { 
        background: transparent !important; 
    }
    .gradio-container input, 
    .gradio-container textarea { 
        background: #1b222c !important; 
        color: #e6edf3 !important; 
        border-color: #2b3440 !important; 
    }
    .gradio-container .btn, 
    .gradio-container button { 
        background: #1b222c !important; 
        color: #e6edf3 !important; 
        border-color: #2b3440 !important; 
    }
    a[aria-label="Use via API"], 
    a[aria-label="Settings"], 
    button[aria-label="Use via API"], 
    button[aria-label="Settings"] { 
        background: #0e1217 !important; 
        color: #e6edf3 !important; 
        border-color: #2b3440 !important; 
    }
    
    /* Footer with reference link */
    #tds-footer { 
        position: fixed; 
        bottom: 10px; 
        left: 12px; 
        z-index: 998; 
        display: flex; 
        align-items: center; 
        gap: 0; 
    }
    #tds-footer a { 
        color: #e6edf3; 
        text-decoration: none; 
        background: #0e1217; 
        padding: 6px 10px; 
        border-radius: 8px; 
        border: 1px solid #2b3440; 
        margin-left: 8px; 
    }
    #tds-footer a:hover { 
        background: #1b222c; 
    }
    #tds-footer svg { 
        width: 16px; 
        height: 16px; 
        vertical-align: middle; 
        margin-right: 6px; 
    }
    </style>
    
    <script>
    // Clear any localStorage/sessionStorage to prevent state persistence issues
    try {
      if (window.localStorage) localStorage.clear();
      if (window.sessionStorage) sessionStorage.clear();
    } catch (e) {}
    
    // Unregister any service workers
    try {
      if (navigator && navigator.serviceWorker) {
        navigator.serviceWorker.getRegistrations().then(rs => {
          rs.forEach(r => { try { r.unregister(); } catch (_) {} });
        }).catch(() => {});
      }
    } catch (_) {}
    
    // Dynamically place reference link near settings button
    (function(){
      const placeRef = () => {
        const settingsEl = document.querySelector('[aria-label="Settings"]');
        const refEl = document.querySelector('#tds-footer a');
        if (settingsEl && refEl && !refEl._moved) {
          settingsEl.parentNode.insertBefore(refEl, settingsEl.nextSibling);
          refEl._moved = true;
          const footer = document.getElementById('tds-footer');
          if (footer) footer.remove();
        }
      };
      const obs = new MutationObserver(placeRef);
      obs.observe(document.body, { childList: true, subtree: true });
      placeRef();
    })();
    </script>
    """)
    
    # ========================================================================
    # Header
    # ========================================================================
    with gr.Row(elem_id="tds-header"):
        gr.Markdown("# TDS: Target Discovery System — Oncology Target Discovery Chatbot")
    
    # ========================================================================
    # Example Queries
    # ========================================================================
    with gr.Accordion("Examples", open=False):
        gr.Markdown("Which target is safer for Lung Cancer: MUC1 or TACSTD2?")
    
    # ========================================================================
    # Search Method Toggle
    # ========================================================================
    with gr.Accordion("Search Method", open=False):
        method_toggle = gr.Radio(
            ["local", "global"], 
            label="Search Method", 
            value="global",
            info="Local: focused neighborhood search. Global: comprehensive graph analysis."
        )
    
    # ========================================================================
    # Error Modal Components
    # ========================================================================
    error_overlay = gr.HTML("", visible=False, elem_id="tds-error-overlay")
    error_md = gr.Markdown("", visible=False, elem_id="tds-error-modal")
    close_btn = gr.Button("×", visible=False, elem_id="tds-error-close")
    
    # ========================================================================
    # Chat Interface
    # ========================================================================
    gr.ChatInterface(
        fn=graphrag_query,
        additional_inputs=[method_toggle],
        additional_outputs=[error_overlay, error_md, close_btn],
        examples=[["Which target is safer for Lung Cancer: MUC1 or TACSTD2?"]],
        title=None,  # Title handled in header
        description=None,
    )
    
    # ========================================================================
    # Footer with Paper Reference
    # ========================================================================
    gr.HTML("""
    <div id="tds-footer">
      <a href="https://www.biorxiv.org/content/10.1101/2025.05.06.652559v1" target="_blank" rel="noopener">
        <svg viewBox="0 0 24 24" fill="#e6edf3" xmlns="http://www.w3.org/2000/svg">
          <path d="M4 3h12l4 4v14a1 1 0 0 1-1 1H4a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1zm12 2H5v16h14V8h-3a2 2 0 0 1-2-2V5zm-7 6h8v2H9v-2zm0 4h8v2H9v-2zm0-8h5v2H9V7z"/>
        </svg>
        Reference Paper
      </a>
    </div>
    """)
    
    # ========================================================================
    # Event Handlers
    # ========================================================================
    close_btn.click(
        lambda: (
            gr.update(visible=False), 
            gr.update(visible=False), 
            gr.update(visible=False)
        ), 
        None, 
        [error_overlay, error_md, close_btn]
    )


# ============================================================================
# Application Entry Point
# ============================================================================

if __name__ == "__main__":
    demo.launch(
        theme=custom_theme,
        server_port=7861,
        server_name="0.0.0.0",  # Accessible from network
        show_error=True,
        quiet=False
    )
