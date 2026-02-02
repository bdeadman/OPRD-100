"""
Unit tests for security validation module.
Tests various attack vectors and edge cases.
"""
import json
import os
import tempfile
import pytest
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from security import SecurityValidator, validate_submission


class TestSecurityValidator:
    """Test suite for SecurityValidator class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.validator = SecurityValidator()
        self.temp_dir = tempfile.mkdtemp()
    
    def create_temp_submission(self, data):
        """Helper to create a temporary submission file."""
        filepath = os.path.join(self.temp_dir, 'test_submission.json')
        with open(filepath, 'w') as f:
            json.dump(data, f)
        return filepath
    
    def test_valid_submission(self):
        """Test that a valid submission passes."""
        valid_submission = {
            "submitter_name": "Test User",
            "reactions": [
                {
                    "Reference": "10.1021/test",
                    "Location": {"Type": "Scheme", "Num": "1"},
                    "Reaction": "C>>CC",
                    "Steps": [
                        {
                            "Reagents": [],
                            "Solvents": [],
                            "Temperature": "25 °C",
                            "Time": "1 h"
                        }
                    ]
                }
            ]
        }
        filepath = self.create_temp_submission(valid_submission)
        is_valid, error, data = validate_submission(filepath)
        assert is_valid, f"Valid submission failed: {error}"
    
    def test_code_injection_eval(self):
        """Test detection of eval() injection."""
        malicious = {
            "submitter_name": "eval(__import__('os').system('ls'))",
            "reactions": []
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Suspicious pattern" in error or "reactions list cannot be empty" in error
    
    def test_code_injection_import(self):
        """Test detection of import statements."""
        malicious = {
            "submitter_name": "Test",
            "reactions": [
                {
                    "Reference": "__import__('subprocess').call(['rm', '-rf', '/'])",
                    "Steps": []
                }
            ]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Suspicious pattern" in error
    
    def test_path_traversal(self):
        """Test detection of path traversal attempts."""
        malicious = {
            "submitter_name": "../../etc/passwd",
            "reactions": [{"Steps": []}]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Suspicious pattern" in error
    
    def test_excessive_nesting(self):
        """Test detection of excessive JSON nesting."""
        # Create deeply nested structure
        nested = {"a": {}}
        current = nested["a"]
        for i in range(15):  # Exceed MAX_NESTING_DEPTH
            current["a"] = {}
            current = current["a"]
        
        malicious = {
            "submitter_name": "Test",
            "reactions": [nested]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "nesting depth" in error.lower()
    
    def test_oversized_string(self):
        """Test detection of oversized strings."""
        malicious = {
            "submitter_name": "A" * 20000,  # Exceed MAX_STRING_LENGTH
            "reactions": [{"Steps": []}]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "String too long" in error or "too long" in error.lower()
    
    def test_too_many_reactions(self):
        """Test detection of excessive reactions."""
        malicious = {
            "submitter_name": "Test",
            "reactions": [{"Steps": []} for _ in range(1500)]  # Exceed MAX_REACTIONS
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Too many reactions" in error
    
    def test_missing_required_keys(self):
        """Test detection of missing required keys."""
        malicious = {
            "submitter_name": "Test"
            # Missing 'reactions' key
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Missing required keys" in error
    
    def test_unexpected_keys(self):
        """Test detection of unexpected keys."""
        malicious = {
            "submitter_name": "Test",
            "reactions": [],
            "malicious_key": "evil_value",
            "backdoor": "payload"
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Unexpected keys" in error
    
    def test_empty_reactions_list(self):
        """Test detection of empty reactions list."""
        malicious = {
            "submitter_name": "Test",
            "reactions": []
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "reactions list cannot be empty" in error
    
    def test_invalid_file_extension(self):
        """Test rejection of non-JSON files."""
        filepath = os.path.join(self.temp_dir, 'test.txt')
        with open(filepath, 'w') as f:
            f.write("not json")
        
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Invalid file extension" in error
    
    def test_nonexistent_file(self):
        """Test handling of nonexistent files."""
        filepath = "/tmp/nonexistent_file_12345.json"
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "does not exist" in error
    
    def test_invalid_json(self):
        """Test handling of malformed JSON."""
        filepath = os.path.join(self.temp_dir, 'invalid.json')
        with open(filepath, 'w') as f:
            f.write("{invalid json content")
        
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Invalid JSON" in error
    
    def test_subprocess_injection(self):
        """Test detection of subprocess attempts."""
        malicious = {
            "submitter_name": "Test",
            "reactions": [
                {
                    "Reference": "subprocess.call(['malicious'])",
                    "Steps": []
                }
            ]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Suspicious pattern" in error
    
    def test_pickle_injection(self):
        """Test detection of pickle deserialization."""
        malicious = {
            "submitter_name": "Test",
            "reactions": [
                {
                    "Reference": "pickle.loads(malicious_data)",
                    "Steps": []
                }
            ]
        }
        filepath = self.create_temp_submission(malicious)
        is_valid, error, _ = validate_submission(filepath)
        assert not is_valid
        assert "Suspicious pattern" in error


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
