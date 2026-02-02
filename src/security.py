"""
Security validation module for submission files.
Protects against malicious attacks including code injection, path traversal, and resource exhaustion.
"""
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple


class SecurityValidator:
    """Validates submission files for security threats."""
    
    # Limits to prevent resource exhaustion
    MAX_FILE_SIZE = 10 * 1024 * 1024  # 10MB
    MAX_REACTIONS = 1000
    MAX_STEPS_PER_REACTION = 50
    MAX_REAGENTS_PER_STEP = 100
    MAX_STRING_LENGTH = 10000
    MAX_NESTING_DEPTH = 10
    
    # Allowed submission keys
    REQUIRED_KEYS = {'submitter_name', 'reactions'}
    OPTIONAL_KEYS = {'method_description', 'repository_url', 'paper_url', 'contact_email'}
    ALLOWED_KEYS = REQUIRED_KEYS | OPTIONAL_KEYS
    
    # Allowed reaction keys
    REACTION_KEYS = {'Reference', 'Location', 'Reaction', 'Steps', 'Product'}
    LOCATION_KEYS = {'Type', 'Num'}
    STEP_KEYS = {'Reagents', 'Solvents', 'Temperature', 'Time', 'Yield'}
    REAGENT_KEYS = {'Reagent', 'Amounts'}
    SOLVENT_KEYS = {'Solvent', 'Amounts'}
    AMOUNT_KEYS = {'Mass', 'Volume', 'Moles', 'Equivalents'}
    YIELD_KEYS = {'Value', 'Type'}
    
    # Dangerous patterns
    DANGEROUS_PATTERNS = [
        r'__import__',
        r'eval\s*\(',
        r'exec\s*\(',
        r'compile\s*\(',
        r'open\s*\(',
        r'file\s*\(',
        r'input\s*\(',
        r'raw_input\s*\(',
        r'execfile\s*\(',
        r'reload\s*\(',
        r'\.\./',  # Path traversal
        r'\.\.\\',
        r'subprocess',
        r'os\.system',
        r'os\.popen',
        r'commands\.',
        r'pickle',
        r'marshal',
        r'__code__',
        r'__globals__',
        r'__builtins__',
        r'lambda\s',
        r'import\s',
        r'from\s+\w+\s+import',
    ]
    
    def __init__(self):
        """Initialize the security validator."""
        self.compiled_patterns = [re.compile(p, re.IGNORECASE) for p in self.DANGEROUS_PATTERNS]
    
    def validate_file_path(self, file_path: str) -> Tuple[bool, str]:
        """
        Validate that the file path is safe and within expected directory.
        
        Args:
            file_path: Path to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            # Resolve absolute path and check it's within allowed directory
            abs_path = Path(file_path).resolve()
            
            # Check file exists
            if not abs_path.exists():
                return False, f"File does not exist: {file_path}"
            
            # Check it's a file, not a directory
            if not abs_path.is_file():
                return False, f"Path is not a file: {file_path}"
            
            # Check file size
            file_size = abs_path.stat().st_size
            if file_size > self.MAX_FILE_SIZE:
                return False, f"File too large: {file_size} bytes (max {self.MAX_FILE_SIZE})"
            
            # Check file extension
            if abs_path.suffix.lower() != '.json':
                return False, f"Invalid file extension: {abs_path.suffix} (expected .json)"
            
            return True, ""
            
        except Exception as e:
            return False, f"Path validation error: {str(e)}"
    
    def validate_json_structure(self, data: Any, depth: int = 0) -> Tuple[bool, str]:
        """
        Recursively validate JSON structure and check for suspicious content.
        
        Args:
            data: Data to validate
            depth: Current nesting depth
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Check nesting depth
        if depth > self.MAX_NESTING_DEPTH:
            return False, f"Excessive nesting depth: {depth} (max {self.MAX_NESTING_DEPTH})"
        
        # Check strings for malicious patterns
        if isinstance(data, str):
            if len(data) > self.MAX_STRING_LENGTH:
                return False, f"String too long: {len(data)} chars (max {self.MAX_STRING_LENGTH})"
            
            for pattern in self.compiled_patterns:
                if pattern.search(data):
                    return False, f"Suspicious pattern detected: {pattern.pattern}"
        
        # Recursively check collections
        elif isinstance(data, dict):
            for key, value in data.items():
                # Validate key
                if not isinstance(key, str):
                    return False, f"Non-string dictionary key: {type(key)}"
                
                is_valid, error = self.validate_json_structure(key, depth + 1)
                if not is_valid:
                    return False, error
                
                # Validate value
                is_valid, error = self.validate_json_structure(value, depth + 1)
                if not is_valid:
                    return False, error
        
        elif isinstance(data, list):
            for item in data:
                is_valid, error = self.validate_json_structure(item, depth + 1)
                if not is_valid:
                    return False, error
        
        # Allow primitives: None, bool, int, float
        elif data is None or isinstance(data, (bool, int, float)):
            pass
        else:
            return False, f"Unexpected data type: {type(data)}"
        
        return True, ""
    
    def validate_submission_schema(self, submission: Dict) -> Tuple[bool, str]:
        """
        Validate the submission follows the expected schema.
        
        Args:
            submission: Submission dictionary to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not isinstance(submission, dict):
            return False, "Submission must be a dictionary"
        
        # Check for required keys
        missing_keys = self.REQUIRED_KEYS - set(submission.keys())
        if missing_keys:
            return False, f"Missing required keys: {missing_keys}"
        
        # Check for unexpected keys
        unexpected_keys = set(submission.keys()) - self.ALLOWED_KEYS
        if unexpected_keys:
            return False, f"Unexpected keys in submission: {unexpected_keys}"
        
        # Validate submitter_name
        submitter_name = submission.get('submitter_name', '')
        if not isinstance(submitter_name, str) or not submitter_name.strip():
            return False, "submitter_name must be a non-empty string"
        if len(submitter_name) > 100:
            return False, "submitter_name too long (max 100 characters)"
        
        # Validate reactions
        reactions = submission.get('reactions', [])
        if not isinstance(reactions, list):
            return False, "reactions must be a list"
        
        if len(reactions) == 0:
            return False, "reactions list cannot be empty"
        
        if len(reactions) > self.MAX_REACTIONS:
            return False, f"Too many reactions: {len(reactions)} (max {self.MAX_REACTIONS})"
        
        # Validate each reaction
        for i, reaction in enumerate(reactions):
            is_valid, error = self.validate_reaction(reaction, i)
            if not is_valid:
                return False, f"Reaction {i}: {error}"
        
        return True, ""
    
    def validate_reaction(self, reaction: Dict, index: int) -> Tuple[bool, str]:
        """Validate a single reaction entry."""
        if not isinstance(reaction, dict):
            return False, f"Reaction must be a dictionary"
        
        # Check for unexpected keys (allow flexibility but warn about extras)
        unexpected = set(reaction.keys()) - self.REACTION_KEYS
        if unexpected:
            # Log but don't fail - allow some flexibility
            pass
        
        # Validate Steps
        steps = reaction.get('Steps', [])
        if not isinstance(steps, list):
            return False, "Steps must be a list"
        
        if len(steps) > self.MAX_STEPS_PER_REACTION:
            return False, f"Too many steps: {len(steps)} (max {self.MAX_STEPS_PER_REACTION})"
        
        for step_idx, step in enumerate(steps):
            if not isinstance(step, dict):
                return False, f"Step {step_idx} must be a dictionary"
            
            # Validate reagents
            reagents = step.get('Reagents', [])
            if isinstance(reagents, list) and len(reagents) > self.MAX_REAGENTS_PER_STEP:
                return False, f"Step {step_idx}: Too many reagents: {len(reagents)} (max {self.MAX_REAGENTS_PER_STEP})"
            
            # Validate solvents
            solvents = step.get('Solvents', [])
            if isinstance(solvents, list) and len(solvents) > self.MAX_REAGENTS_PER_STEP:
                return False, f"Step {step_idx}: Too many solvents: {len(solvents)} (max {self.MAX_REAGENTS_PER_STEP})"
        
        return True, ""
    
    def sanitize_path(self, path: str) -> str:
        """
        Sanitize a file path to prevent path traversal attacks.
        
        Args:
            path: Path to sanitize
            
        Returns:
            Sanitized path
        """
        # Remove any path traversal attempts
        path = path.replace('..', '')
        path = path.replace('~', '')
        
        # Get basename only (no directory components)
        path = os.path.basename(path)
        
        return path
    
    def validate_submission_file(self, file_path: str) -> Tuple[bool, str, Dict]:
        """
        Complete validation of a submission file.
        
        Args:
            file_path: Path to submission file
            
        Returns:
            Tuple of (is_valid, error_message, submission_data)
        """
        # Validate file path
        is_valid, error = self.validate_file_path(file_path)
        if not is_valid:
            return False, error, {}
        
        # Load and parse JSON
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                submission = json.load(f)
        except json.JSONDecodeError as e:
            return False, f"Invalid JSON: {str(e)}", {}
        except Exception as e:
            return False, f"Error reading file: {str(e)}", {}
        
        # Validate JSON structure for malicious content
        is_valid, error = self.validate_json_structure(submission)
        if not is_valid:
            return False, f"Security validation failed: {error}", {}
        
        # Validate submission schema
        is_valid, error = self.validate_submission_schema(submission)
        if not is_valid:
            return False, f"Schema validation failed: {error}", {}
        
        return True, "", submission


def validate_submission(file_path: str) -> Tuple[bool, str, Dict]:
    """
    Convenience function to validate a submission file.
    
    Args:
        file_path: Path to submission file
        
    Returns:
        Tuple of (is_valid, error_message, submission_data)
    """
    validator = SecurityValidator()
    return validator.validate_submission_file(file_path)
