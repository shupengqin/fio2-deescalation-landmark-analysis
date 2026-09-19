"""Data-free syntax and publication-boundary checks; does not run analyses."""
from pathlib import Path
import ast,json,re
ROOT=Path(__file__).resolve().parents[1]
allowed={'.py','.sql','.md','.txt','.json'}
errors=[];count=0
for p in ROOT.rglob('*'):
    if not p.is_file() or '.git' in p.parts or '__pycache__' in p.parts:
        continue
    rel=p.relative_to(ROOT)
    if p.name!='.gitignore' and p.suffix not in allowed:
        errors.append(f'Unexpected release file: {rel}')
        continue
    text=p.read_text(encoding='utf8')
    if p.suffix=='.py':
        ast.parse(text,filename=str(rel));count+=1
    for pattern in [r'gh[pousr]_[A-Za-z0-9]{20,}', r'github_pat_[A-Za-z0-9_]{20,}',r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----']:
        if re.search(pattern,text):errors.append(f'Potential credential: {rel}')
    if p.suffix in {'.py','.sql'} and p.name!='validate_release.py':
        if re.search(r"(?:password|token|api_key)\s*=\s*['\"][^'\"]+['\"]",text,re.I):errors.append(f'Literal credential assignment: {rel}')
        if re.search(r'[A-Z]:\\(?:Users|SICDB)\\',text):errors.append(f'Unportable source path: {rel}')
assert not errors,'\n'.join(errors)
print(json.dumps({'python_files_parsed':count,'release_file_boundary':'passed','credential_pattern_scan':'passed','patient_data_reanalysis':False}))
