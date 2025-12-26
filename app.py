import streamlit as st
import json
import pandas as pd
import os

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="ReaxFF Library",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- LOADING DATA ---
@st.cache_data
def load_database():
    file_path = "data/master_index.json"
    if not os.path.exists(file_path):
        return []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error loading database: {e}")
        return []

def get_snippet(filename, elements):
    """Generates a LAMMPS input snippet."""
    elems_str = " ".join(elements)
    return f"""# LAMMPS Input Snippet
pair_style reax/c lmp_control
pair_coeff * * {filename} {elems_str}
fix qeq all qeq/reax 1 0.0 10.0 1e-6 param.qeq"""

# --- MAIN APP LOGIC ---
def main():
    st.title("⚛️ ReaxFF Potential Library")
    st.markdown("""
    **Interactive database for Reactive Force Fields.** Scraped from public repositories (GitHub) and curated automatically.
    """)

    # Load Data
    raw_data = load_database()
    
    if not raw_data:
        st.warning("Database is empty. Please run the crawler and cleaner scripts first.")
        st.stop()

    df = pd.DataFrame(raw_data)

    # --- SIDEBAR: FILTERS ---
    st.sidebar.header("🔍 Search Filters")
    
    # 1. Element Selector
    # Get unique elements from the entire dataset
    all_elements = set()
    for el_list in df['elements']:
        all_elements.update(el_list)
    sorted_elements = sorted(list(all_elements))

    selected_elements = st.sidebar.multiselect(
        "Required Elements:",
        options=sorted_elements,
        default=["C", "H", "O"],
        help="Select the elements that MUST be present in the potential."
    )

    # 2. Strict Mode
    strict_mode = st.sidebar.checkbox(
        "Strict Match", 
        value=False,
        help="If checked, finds potentials containing EXACTLY the selected elements (no more, no less)."
    )

    # --- FILTERING LOGIC ---
    if selected_elements:
        req_set = set(selected_elements)
        
        # Apply filter
        def filter_func(row_elements):
            row_set = set(row_elements)
            if strict_mode:
                return row_set == req_set
            else:
                return req_set.issubset(row_set)

        filtered_df = df[df['elements'].apply(filter_func)]
    else:
        filtered_df = df

    # --- RESULTS DISPLAY ---
    st.divider()
    col_info, col_count = st.columns([3, 1])
    col_info.subheader("Available Potentials")
    col_count.metric("Count", len(filtered_df))

    if not filtered_df.empty:
        # Display Table
        st.dataframe(
            filtered_df[['system', 'original_filename', 'source_repo']],
            use_container_width=True,
            column_config={
                "system": "System (Elements)",
                "original_filename": "Filename",
                "source_repo": "Source Repository"
            },
            hide_index=True,
            selection_mode="single-row",
            on_select="rerun" # Allows interactive selection (Streamlit newer versions)
        )

        # --- SELECTION & DETAILS ---
        # User selects a file to download/view
        st.subheader("📂 File Details & Download")
        
        selected_file_id = st.selectbox(
            "Select a file to inspect:",
            options=filtered_df['id'],
            format_func=lambda x: filtered_df[filtered_df['id'] == x]['original_filename'].values[0]
        )

        if selected_file_id:
            record = filtered_df[filtered_df['id'] == selected_file_id].iloc[0]
            
            c1, c2 = st.columns(2)
            
            with c1:
                st.info(f"**Source:** {record['source_repo']}")
                st.text(f"Path: {record['local_path']}")
                if record.get('download_url'):
                    st.link_button("View on GitHub", record['download_url'])

            with c2:
                # Code Snippet
                snippet = get_snippet(record['original_filename'], record['elements'])
                st.code(snippet, language="bash")

            # Content Preview and Download
            # Try to read local file
            local_path = record['local_path']
            # Fix path for Windows/Linux compatibility if running on cloud
            # Cloud runs on Linux, your path might have backslashes. Let's fix it.
            clean_path = local_path.replace("\\", "/") 
            
            if os.path.exists(clean_path):
                with open(clean_path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                with st.expander("📄 View File Content (Preview)"):
                    st.text(content[:2000] + "\n... (truncated)")

                st.download_button(
                    label=f"⬇️ Download {record['original_filename']}",
                    data=content,
                    file_name=record['original_filename'],
                    mime="text/plain",
                    type="primary"
                )
            else:
                st.error(f"File not found on server at: {clean_path}")

    else:
        st.info("No potentials found matching these criteria. Try removing some elements.")

if __name__ == "__main__":
    main()