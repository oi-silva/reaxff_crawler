from abc import ABC, abstractmethod
from typing import List, Dict, Any

class DataSourceInterface(ABC):
    """Base contract for data connectors."""

    @abstractmethod
    def authenticate(self) -> None:
        pass

    @abstractmethod
    def search(self, query_term: str) -> List[Dict[str, Any]]:
        """Returns file metadata list: [{'name': '...', 'download_url': '...'}]"""
        pass

    @abstractmethod
    def get_file_content(self, download_url: str) -> str:
        pass