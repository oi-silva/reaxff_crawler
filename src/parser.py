import hashlib
import re
from typing import List, Dict, Optional, Any

class ReaxFFParser:
    """
    Responsible for validating and extracting chemical metadata from ReaxFF files.
    """

    @staticmethod
    def get_content_hash(content: str) -> str:
        return hashlib.md5(content.encode('utf-8')).hexdigest()

    @staticmethod
    def parse(content: str) -> Dict[str, Any]:
        lines = content.splitlines()
        result = {
            "valid": False,
            "atoms": [],
            "error": None
        }

        if not lines:
            result["error"] = "Empty file"
            return result

        try:
            # --- STRATEGY: Flexible Search ---
            is_reax = False
            atom_count = 0
            atom_line_index = -1

            # We scan the first 300 lines looking for the atom count line.
            # Pattern: Starts with an Integer, followed by whitespace, then '!' (optional), 
            # and containing the word "atom" (case insensitive).
            
            for i, line in enumerate(lines[:300]):
                clean_line = line.strip().lower()
                
                # Skip comments or empty lines at the very top
                if not clean_line or clean_line.startswith("!"):
                    continue

                parts = clean_line.split()
                
                # Check 1: Does it start with a number?
                if parts and parts[0].isdigit():
                    # Check 2: Does the rest of the line contain "atom"?
                    rest_of_line = " ".join(parts[1:])
                    if "atom" in rest_of_line:
                        atom_count = int(parts[0])
                        atom_line_index = i
                        is_reax = True
                        break

            if not is_reax:
                # Debugging aid: print first 5 lines to see what's wrong if needed
                # print(f"DEBUG HEAD: {lines[:5]}")
                result["error"] = "Could not find a line matching pattern: [INT] ... 'atom' ..."
                return result

            # --- Extract Elements ---
            atoms = []
            current_line = atom_line_index + 1
            
            # Read 'atom_count' lines immediately following the declaration
            while len(atoms) < atom_count and current_line < len(lines):
                line = lines[current_line].strip()
                
                # Skip pure comments or empty lines
                if not line or line.startswith("!"):
                    current_line += 1
                    continue
                
                parts = line.split()
                if parts:
                    # Valid ReaxFF element lines start with the Element Symbol (e.g., "C  12.00")
                    element = parts[0]
                    # Simple check: Elements are usually 1 or 2 letters. 
                    # Excluding numbers prevents reading parameters as elements.
                    if element.isalpha() and len(element) <= 2:
                        atoms.append(element)
                    elif len(atoms) > 0:
                        # If we hit a number after reading some atoms, we might be done
                        # (Some files don't have exactly atom_count lines if structured weirdly)
                        pass
                
                current_line += 1

            # Validation logic
            if len(atoms) > 0: 
                # We relax the strict count check because sometimes duplicate lines exist or comments intervene
                result["valid"] = True
                result["atoms"] = sorted(list(set(atoms))) # Deduplicate and sort
            else:
                 result["error"] = f"Found atom count {atom_count} but could not extract element symbols."

        except Exception as e:
            result["valid"] = False
            result["error"] = f"Parsing exception: {str(e)}"

        return result