import os
from dotenv import load_dotenv
from src.connectors.github_connector import GitHubConnector
from src.file_manager import FileManager

load_dotenv()

# --- CONFIGURATION ---
# Set to True to download EVERYTHING found. Set to False to limit for testing.
FULL_SCRAP_MODE = True 

# Limits (Only used if FULL_SCRAP_MODE is False)
TEST_SEARCH_PAGES = 1      # How many pages of search results to fetch
TEST_DOWNLOAD_LIMIT = 5    # How many files to actually download

def main():
    print("--- Starting ReaxFF Crawler ---")
    print(f"[CONFIG] Mode: {'FULL SCRAP' if FULL_SCRAP_MODE else 'TESTING'}")
    
    # 1. Init
    token = os.getenv("GITHUB_TOKEN")
    github_source = GitHubConnector(api_token=token)
    file_manager = FileManager(base_path="data/raw")
    
    # 2. Search
    query = "ffield.reax"
    # Logic: If Full Mode, fetch 5 pages (500 items). If Test, fetch 1 page.
    # Note: GitHub limits code search results to ~1000 items total anyway.
    pages_to_fetch = 10 if FULL_SCRAP_MODE else TEST_SEARCH_PAGES
    
    results = github_source.search(query, max_pages=pages_to_fetch)
    
    # 3. Download Pipeline
    download_count = 0
    # Logic: If Full Mode, download all results. If Test, stop at limit.
    limit = len(results) if FULL_SCRAP_MODE else TEST_DOWNLOAD_LIMIT
    
    print(f"\n[ACTION] Starting download pipeline (Target: {limit} files)...")

    for i, item in enumerate(results):
        if i >= limit:
            print("[INFO] Test limit reached. Stopping downloads.")
            break

        # Check existing (Using dictionary fields)
        if file_manager.file_exists(item['source'], item['repo'], item['name']):
            print(f" -> Skipping ({i+1}/{limit}): {item['name']} (Already exists)")
            continue

        print(f" -> Downloading ({i+1}/{limit}): {item['name']}")
        content = github_source.get_file_content(item['download_url'])
        
        if content:
            # CHANGED: Now passing the full 'item' dictionary as metadata
            saved_path = file_manager.save_file(
                content=content,
                metadata=item
            )
            
            if saved_path:
                print(f"    [SUCCESS] Saved.")
                download_count += 1
            else:
                print(f"    [FAILURE] Save error.")
        else:
            print(f"    [FAILURE] Empty content.")

    print(f"\n--- Process Complete. Downloaded {download_count} new files. ---")

if __name__ == "__main__":
    main()