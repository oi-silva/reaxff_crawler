import requests
import time
from typing import List, Dict, Any
from src.interfaces import DataSourceInterface

class GitHubConnector(DataSourceInterface):
    def __init__(self, api_token: str = None):
        self.base_url = "https://api.github.com"
        self.api_token = api_token
        self.headers = {"Accept": "application/vnd.github.v3+json"}
        self.authenticate()

    def authenticate(self) -> None:
        if self.api_token:
            self.headers["Authorization"] = f"token {self.api_token}"
        else:
            print("[WARNING] No token provided. Rate limits will be strict.")

    def search(self, query_term: str, max_pages: int = 1) -> List[Dict[str, Any]]:
        """
        Searches GitHub with pagination support.
        Args:
            query_term: The search keyword.
            max_pages: Limit of pages to fetch (1 page = ~100 items). Set higher for full scraping.
        """
        results = []
        page = 1
        per_page = 100  # GitHub API maximum per page
        
        print(f"[INFO] Starting search for '{query_term}' (Max pages: {max_pages})...")

        while page <= max_pages:
            search_url = f"{self.base_url}/search/code?q=filename:{query_term}&per_page={per_page}&page={page}"
            
            try:
                # Sleep to respect rate limits (essential for bulk scraping)
                if page > 1:
                    time.sleep(2.0)

                response = requests.get(search_url, headers=self.headers)
                
                if response.status_code == 200:
                    data = response.json()
                    items = data.get("items", [])
                    
                    if not items:
                        print(f"[INFO] Page {page}: No more items found.")
                        break

                    print(f"[INFO] Page {page}: Found {len(items)} items.")

                    for item in items:
                        meta = {
                            "name": item.get("name"),
                            "path": item.get("path"),
                            "repo": item.get("repository", {}).get("full_name"),
                            "download_url": item.get("html_url").replace("github.com", "raw.githubusercontent.com").replace("/blob/", "/"),
                            "source": "GitHub"
                        }
                        results.append(meta)
                    
                    page += 1
                    
                elif response.status_code == 403:
                    print("[WARNING] Rate limit hit (403). Stopping search early.")
                    break
                else:
                    print(f"[ERROR] API Status {response.status_code}: {response.text}")
                    break

            except Exception as e:
                print(f"[ERROR] Search exception on page {page}: {e}")
                break

        print(f"[INFO] Search complete. Total candidates found: {len(results)}")
        return results

    def get_file_content(self, download_url: str) -> str:
        try:
            time.sleep(0.5) 
            response = requests.get(download_url)
            if response.status_code == 200:
                return response.text
        except Exception as e:
            print(f"[ERROR] Download failed: {e}")
        return ""