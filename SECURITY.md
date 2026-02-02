# Security Measures

This document outlines the security measures implemented to protect the OPRD-100 repository from malicious submissions.

## Overview

The repository implements multiple layers of security to prevent:
- **Code injection attacks** - Malicious code in JSON submissions
- **Path traversal attacks** - Accessing files outside submission directory
- **Resource exhaustion attacks** - Large files or excessive computation
- **Command injection** - Shell commands in submission data
- **Denial of Service** - Runaway processes or infinite loops

## Security Layers

### 1. File Path Validation (GitHub Actions)
- Only processes files in `data/submissions/` directory
- Validates file extension is `.json`
- Checks file size (max 10MB)
- Sanitizes submitter names (alphanumeric only)

### 2. Content Validation (Python)
The `security.py` module validates all submission content:

**File Limits:**
- Maximum file size: 10MB
- Maximum reactions: 1,000
- Maximum steps per reaction: 50
- Maximum reagents per step: 100
- Maximum string length: 10,000 characters
- Maximum JSON nesting depth: 10 levels

**Pattern Detection:**
Scans for dangerous patterns including:
- Python code execution (`eval`, `exec`, `__import__`)
- File operations (`open`, `file`)
- System commands (`subprocess`, `os.system`)
- Path traversal (`../`, `..\`)
- Serialization exploits (`pickle`, `marshal`)
- Lambda functions and imports

**Schema Validation:**
- Enforces required fields
- Rejects unexpected keys
- Validates data types
- Checks for empty or malformed data

### 3. Resource Limits (GitHub Actions)
- **Job timeout**: 30 minutes maximum
- **Virtual memory**: 4GB limit
- **CPU time**: 30 minutes limit
- **Process timeout**: 25 minutes with forced termination

### 4. JSON Security
- Only primitive types allowed (string, number, boolean, null, dict, list)
- No code execution during JSON parsing
- Recursive validation of nested structures
- Type checking for all values

### 5. Input Sanitization
All user-provided strings are sanitized:
- Submitter names limited to alphanumeric + spaces/hyphens/underscores
- Maximum length enforcement
- Pattern matching against dangerous content
- Path components stripped from filenames

## Usage

### Validating a Submission File

```python
from src.security import validate_submission

# Validate submission
is_valid, error_msg, submission_data = validate_submission('path/to/submission.json')

if not is_valid:
    print(f"Validation failed: {error_msg}")
    exit(1)

# Proceed with validated data
print(f"Submission valid: {submission_data['submitter_name']}")
```

### Running Scoring (Automatic Validation)

The `run_scoring.py` script automatically validates submissions:

```bash
python scripts/run_scoring.py \
    --submission-file data/submissions/my_submission.json \
    --output-dir results/my_output
```

If validation fails, the script exits with code 1 and displays the error.

## Security Best Practices

### For Maintainers
1. **Review all PRs manually** before merging submissions
2. **Monitor GitHub Actions logs** for suspicious activity
3. **Keep dependencies updated** to patch security vulnerabilities
4. **Limit repository permissions** to trusted contributors only
5. **Enable branch protection** on main branch

### For Contributors
1. **Test submissions locally** before creating a PR
2. **Follow the submission schema** exactly
3. **Keep file sizes reasonable** (< 1MB recommended)
4. **Report security concerns** via private disclosure

## Incident Response

If a security issue is discovered:

1. **Do NOT merge** the malicious PR
2. **Close the PR** with a comment explaining the issue
3. **Report to GitHub** if necessary
4. **Review logs** for any successful exploits
5. **Update security measures** to prevent future attacks
6. **Document the incident** for future reference

## Security Boundaries

### What IS Protected
✅ Code injection via JSON content  
✅ Path traversal to access system files  
✅ Resource exhaustion (CPU, memory, disk)  
✅ Command injection via submission data  
✅ Malformed JSON or excessive nesting  

### What IS NOT Protected
❌ Network-based attacks (DDoS to GitHub)  
❌ Social engineering attacks  
❌ Compromised GitHub accounts  
❌ Vulnerabilities in dependencies (RDKit, pandas, etc.)  

## Testing Security

Run security validation tests:

```bash
python -m pytest tests/test_security.py -v
```

Test with a malicious submission:

```bash
python scripts/run_scoring.py \
    --submission-file tests/fixtures/malicious_submission.json \
    --output-dir /tmp/test_output
# Should fail with security error
```

## Reporting Security Issues

If you discover a security vulnerability:

1. **Do NOT open a public issue**
2. **Email the maintainers** directly
3. **Provide details** about the vulnerability
4. **Wait for confirmation** before public disclosure

## Security Updates

- **2026-02-02**: Initial security implementation
  - Added SecurityValidator class
  - Integrated validation into scoring pipeline
  - Added GitHub Actions security checks
  - Implemented resource limits and timeouts

## References

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [GitHub Security Best Practices](https://docs.github.com/en/actions/security-guides/security-hardening-for-github-actions)
- [Python Security Best Practices](https://python.readthedocs.io/en/stable/library/security_warnings.html)
