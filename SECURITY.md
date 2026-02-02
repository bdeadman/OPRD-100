# Security Measures

This document outlines the security measures implemented to protect the OPRD-100 repository from malicious submissions.

## Overview

The repository implements multiple layers of security to prevent:
- **Code injection attacks** - Malicious code in JSON submissions
- **Path traversal attacks** - Accessing files outside submission directory
- **Resource exhaustion attacks** - Large files or excessive computation
- **Command injection** - Shell commands in submission data
- **Denial of Service** - Runaway processes or infinite loops
- **Unauthorized file modifications** - PRs that modify code, workflows, or data

## Security Layers

### 1. **PR File Restrictions** (NEW - GitHub Actions)
**Enforces that PRs can ONLY modify files in `data/submissions/`**

A dedicated validation workflow (`validate_pr_changes.yml`) runs on every PR:
- ✅ Scans all changed files in the PR
- ✅ Rejects PRs that modify ANY files outside `data/submissions/`
- ✅ Posts clear error messages explaining what files are not allowed
- ✅ Blocks scoring workflow from running until validation passes
- ✅ Provides helpful guidance on fixing the PR

**Example blocked files:**
- `.github/workflows/*` - Workflow files
- `scripts/*` - Python scripts
- `src/*` - Source code
- `README.md` - Documentation
- `requirements.txt` - Dependencies
- Any other repository files

**This prevents attackers from:**
- Modifying workflows to bypass security
- Injecting malicious code into scoring scripts
- Tampering with validation logic
- Changing dependencies to malicious packages
- Modifying documentation or data

### 2. File Path Validation (GitHub Actions)
- Only processes files in `data/submissions/` directory
- Validates file extension is `.json`
- Checks file size (max 10MB)
- Sanitizes submitter names (alphanumeric only)

### 3. Content Validation (Python)
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

### 4. Resource Limits (GitHub Actions)
- **Job timeout**: 30 minutes maximum
- **Virtual memory**: 4GB limit
- **CPU time**: 30 minutes limit
- **Process timeout**: 25 minutes with forced termination

### 5. JSON Security
- Only primitive types allowed (string, number, boolean, null, dict, list)
- No code execution during JSON parsing
- Recursive validation of nested structures
- Type checking for all values

### 6. Input Sanitization
All user-provided strings are sanitized:
- Submitter names limited to alphanumeric + spaces/hyphens/underscores
- Maximum length enforcement
- Pattern matching against dangerous content
- Path components stripped from filenames

### 7. CODEOWNERS Protection
GitHub CODEOWNERS file requires repository owner approval for:
- Workflow files (`.github/workflows/*`)
- Python scripts (`scripts/*`)
- Source code (`src/*`)
- Dependencies (`requirements*.txt`)
- Security documentation (`SECURITY.md`)
- README and other critical files

## How It Works

### Pull Request Flow

```
1. User creates PR with files
   ↓
2. validate_pr_changes.yml runs
   ↓ (if files outside data/submissions/)
3. ❌ PR BLOCKED - Error comment posted
   User must fix PR
   
   ↓ (if only data/submissions/ files)
4. ✅ Validation passes
   ↓
5. score_submission.yml runs
   ↓
6. Double-checks file restrictions
   ↓
7. Security validates JSON content
   ↓
8. Scoring proceeds (if all checks pass)
```

### What Contributors Can Do
✅ Add new files to `data/submissions/`  
✅ Modify their own submission files in `data/submissions/`  
✅ Delete their own submission files

### What Contributors CANNOT Do
❌ Modify workflow files  
❌ Change Python scripts  
❌ Edit source code  
❌ Update dependencies  
❌ Modify documentation  
❌ Access files outside `data/submissions/`  

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
5. **Enable branch protection** on main branch with required status checks
6. **Update CODEOWNERS** file with actual GitHub usernames
7. **Require "Validate PR File Changes" workflow to pass** before merging

### For Contributors
1. **Only modify files in `data/submissions/`** - anything else will be rejected
2. **Test submissions locally** before creating a PR
3. **Follow the submission schema** exactly
4. **Keep file sizes reasonable** (< 1MB recommended)
5. **Report security concerns** via private disclosure
6. **Do not attempt to bypass restrictions** - violations may result in being blocked

## Setting Up Repository Protection

To fully enable these security measures:

1. **Enable Branch Protection** on `main`:
   - Go to Settings → Branches → Add rule
   - Branch name pattern: `main`
   - Enable: "Require status checks to pass before merging"
   - Select: "Validate PR File Changes" as required check
   - Enable: "Require review from Code Owners"
   - Enable: "Do not allow bypassing the above settings"

2. **Update CODEOWNERS**:
   - Replace `@OWNER_USERNAME` with your GitHub username
   - Add additional maintainers if needed

3. **Set Repository Permissions**:
   - Settings → General → Permissions
   - Disable: "Allow merge commits" (optional, forces clean history)
   - Disable: "Allow forking" (optional, for private repos)

## Incident Response

If a security issue is discovered:

1. **Do NOT merge** the malicious PR
2. **Close and lock the PR** with a comment explaining the issue
3. **Block the user** if intentional malicious activity is detected
4. **Report to GitHub** if necessary
5. **Review logs** for any successful exploits
6. **Update security measures** to prevent future attacks
7. **Document the incident** for future reference

## Security Boundaries

### What IS Protected
✅ Code injection via JSON content  
✅ Path traversal to access system files  
✅ Resource exhaustion (CPU, memory, disk)  
✅ Command injection via submission data  
✅ Malformed JSON or excessive nesting  
✅ **Unauthorized file modifications (NEW)**  
✅ **Workflow tampering (NEW)**  
✅ **Code modification attempts (NEW)**  

### What IS NOT Protected
❌ Network-based attacks (DDoS to GitHub)  
❌ Social engineering attacks  
❌ Compromised GitHub accounts with write access  
❌ Vulnerabilities in dependencies (RDKit, pandas, etc.)  
❌ Zero-day exploits in GitHub Actions  

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

Test file restriction workflow locally:
```bash
# Try to modify a protected file
echo "malicious" >> scripts/run_scoring.py
git add scripts/run_scoring.py
git commit -m "Test"
# This would be blocked by the validation workflow in a PR
```

## Reporting Security Issues

If you discover a security vulnerability:

1. **Do NOT open a public issue**
2. **Email the maintainers** directly
3. **Provide details** about the vulnerability
4. **Wait for confirmation** before public disclosure

## Security Updates

- **2026-02-02**: File restriction enforcement added
  - Added PR file validation workflow
  - Added double-check in scoring workflow
  - Created CODEOWNERS file
  - Updated documentation
  
- **2026-02-02**: Initial security implementation
  - Added SecurityValidator class
  - Integrated validation into scoring pipeline
  - Added GitHub Actions security checks
  - Implemented resource limits and timeouts

## References

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [GitHub Security Best Practices](https://docs.github.com/en/actions/security-guides/security-hardening-for-github-actions)
- [Python Security Best Practices](https://python.readthedocs.io/en/stable/library/security_warnings.html)
- [GitHub Branch Protection](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches)
- [GitHub CODEOWNERS](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-code-owners)
