import os
import json
from pathlib import Path
from src.parser import ReaxFFParser

def main():
    print("--- Starting Data Cleaning & Indexing ---")

    raw_dir = Path("data/raw")
    output_file = Path("data/master_index.json")
    
    unique_potentials = {} # Key: Hash, Value: Metadata
    stats = {"processed": 0, "valid": 0, "duplicates": 0, "invalid": 0}

    # 1. Traverse all folders in raw data
    # We walk through the directory tree
    for root, dirs, files in os.walk(raw_dir):
        for filename in files:
            # We only process the potential files, not the sidecar .json files
            if filename.endswith(".json"):
                continue

            stats["processed"] += 1
            file_path = Path(root) / filename
            
            # Load sidecar metadata if available
            sidecar_path = file_path.with_suffix(file_path.suffix + ".json")
            if not sidecar_path.exists():
                # Fallback: try appending .json to full name (e.g. file.reax -> file.reax.json)
                sidecar_path = Path(str(file_path) + ".json")
            
            meta_origin = {}
            if sidecar_path.exists():
                try:
                    with open(sidecar_path, 'r', encoding='utf-8') as f:
                        meta_origin = json.load(f)
                except:
                    pass # Ignore sidecar errors

            # 2. Read and Parse content
            try:
                with open(file_path, "r", encoding="utf-8", errors='ignore') as f:
                    content = f.read()
            except Exception as e:
                print(f"[ERROR] Could not read {filename}: {e}")
                continue

            # Generate Hash (Deduplication ID)
            file_hash = ReaxFFParser.get_content_hash(content)

            # Check duplication
            if file_hash in unique_potentials:
                stats["duplicates"] += 1
                # Optional: We could track multiple sources for the same file here
                continue

            # Parse Physics/Chemistry
            parse_result = ReaxFFParser.parse(content)

            if parse_result["valid"]:
                stats["valid"] += 1
                
                # Create the clean entry
                entry = {
                    "id": file_hash,
                    "elements": parse_result["atoms"],
                    "system": "-".join(parse_result["atoms"]), # e.g. "C-H-O"
                    "original_filename": filename,
                    "source_repo": meta_origin.get("repo", "unknown"),
                    "local_path": str(file_path),
                    "download_url": meta_origin.get("download_url", "")
                }
                
                unique_potentials[file_hash] = entry
                print(f"[VALID] {entry['system']:<15} | {filename}")
            else:
                stats["invalid"] += 1
                print(f"[INVALID] {filename} - {parse_result['error']}")

    # 3. Save Master Index
    print("\n--- Summary ---")
    print(f"Total Files Scanned: {stats['processed']}")
    print(f"Valid ReaxFF:        {stats['valid']}")
    print(f"Duplicates Ignored:  {stats['duplicates']}")
    print(f"Invalid/Trash:       {stats['invalid']}")
    
    with open(output_file, "w", encoding="utf-8") as f:
        # Convert dict to list for easier JSON consumption
        json.dump(list(unique_potentials.values()), f, indent=4)
    
    print(f"\n[SUCCESS] Master index saved to: {output_file}")

if __name__ == "__main__":
    main()