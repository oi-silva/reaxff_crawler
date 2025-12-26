import json
import sys
from pathlib import Path
from typing import List, Dict

def load_db(db_path: str) -> List[Dict]:
    path = Path(db_path)
    if not path.exists():
        print(f"[ERROR] Database file not found at: {path}")
        return []
    
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)

def search_potentials(required_elements: List[str], database: List[Dict]) -> List[Dict]:
    """
    Finds potentials that contain ALL the required elements.
    Logic: Required_Set must be a SUBSET of Potential_Elements_Set.
    """
    matches = []
    req_set = set(e.upper() for e in required_elements)

    for entry in database:
        available_elements = set(e.upper() for e in entry['elements'])
        
        # Check if the potential has everything we need
        if req_set.issubset(available_elements):
            matches.append(entry)
            
    return matches

def main():
    print("--- ReaxFF Bank Query Tool ---")
    
    # 1. Load Database
    db_file = "data/master_index.json"
    data = load_db(db_file)
    print(f"[INFO] Loaded {len(data)} potentials from index.")

    # 2. Get User Input
    if len(sys.argv) > 1:
        # Input from command line arguments
        query_input = sys.argv[1:]
    else:
        # Input from interactive prompt
        print("\nEnter chemical elements separated by space (e.g., C H O):")
        user_input = input("> ").strip()
        if not user_input:
            print("[INFO] No input provided. Exiting.")
            return
        query_input = user_input.split()

    # Clean input
    required = [x.strip() for x in query_input if x.strip()]
    print(f"\n[ACTION] Searching for potentials containing: {', '.join(required)}...")

    # 3. Perform Search
    results = search_potentials(required, data)

    # 4. Display Results
    if results:
        print(f"\n[SUCCESS] Found {len(results)} matching potentials:\n")
        print(f"{'ID (Hash)':<10} | {'System':<20} | {'Source Repo':<30} | {'Filename'}")
        print("-" * 100)
        
        for res in results:
            short_hash = res['id'][:8]
            system_str = res['system'][:20]
            repo_str = res['source_repo'][:30]
            print(f"{short_hash:<10} | {system_str:<20} | {repo_str:<30} | {res['original_filename']}")
            
        print("\n[TIP] Check 'data/master_index.json' for full paths and URLs.")
    else:
        print("\n[RESULT] No potentials found containing ALL those elements.")
        print("Try removing some elements to broaden your search.")

if __name__ == "__main__":
    main()