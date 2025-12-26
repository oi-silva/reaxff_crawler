import os
import json
from pathlib import Path
from typing import Dict, Any

class FileManager:
    """Handles file system operations."""

    def __init__(self, base_path: str = "data/raw"):
        self.base_path = Path(base_path)
        self._ensure_base_directory()

    def _ensure_base_directory(self) -> None:
        self.base_path.mkdir(parents=True, exist_ok=True)

    def file_exists(self, source: str, repo_name: str, filename: str) -> bool:
        safe_repo = repo_name.replace("/", "_").replace("\\", "_")
        target_file = self.base_path / source / safe_repo / filename
        return target_file.exists()

    def save_file(self, content: str, metadata: Dict[str, Any]) -> str:
        """
        Saves the raw content AND a sidecar JSON metadata file.
        """
        source = metadata.get('source', 'Unknown')
        repo_name = metadata.get('repo', 'Unknown')
        filename = metadata.get('name', 'unknown_file')

        # 1. Path construction
        safe_repo = repo_name.replace("/", "_").replace("\\", "_")
        target_dir = self.base_path / source / safe_repo
        target_dir.mkdir(parents=True, exist_ok=True)

        file_path = target_dir / filename
        json_path = target_dir / f"{filename}.json"

        try:
            # 2. Save Raw Content
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(content)

            # 3. Save Metadata (Sidecar)
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=4)

            return str(file_path)
        except Exception as e:
            print(f"[ERROR] Could not save files for {filename}: {e}")
            return ""