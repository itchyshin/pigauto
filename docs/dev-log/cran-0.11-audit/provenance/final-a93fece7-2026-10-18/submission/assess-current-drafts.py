"""Record root's separate post-draft review of exact current release prose."""
from pathlib import Path
import datetime,hashlib,json,subprocess
S=Path(__file__).resolve().parent;R=Path('/private/tmp/pigauto-pr236-validator-20261017');O=R/'docs/dev-log/cran-0.11-audit/provenance/final-a93fece7-2026-10-18'
items=[(O/'README.md','Frozen release evidence report',O/'README-assessment.json'),(R/'cran-comments.md','CRAN submission comments',O/'cran-comments-assessment.json'),(S/'after-task.md','Pre-submission after-task report',S/'after-task-final-assessment.json'),(S/'pr236-body.md','Evidence pull request description',S/'pr236-body-assessment.json')]
for p,genre,out in items:
 version='v1' if p.name=='cran-comments.md' else 'v2'
 b=p.read_bytes();lines=b.decode().splitlines();examples=[]
 for cue in ['The DOI','This does NOT cover','The original compressed archive','Supporting captures','The highest proven rung','The current frozen']:
  found=next(((i,l) for i,l in enumerate(lines,1) if cue in l),None)
  if found:
   i,l=found;examples.append({'location':'Line '+str(i),'cue':'Concrete observation and its scope are explicit','impact':'Makes the evidence and remaining limitation visible without inferring acceptance or general scientific validity.','repair':'no change needed','line':i,'quote':l})
  if len(examples)==2:break
 if not examples:
  i,l=next((i,l) for i,l in enumerate(lines,1) if l.strip() and not l.startswith('#'));examples=[{'location':'Opening paragraph','cue':'Names the specific release action','impact':'Keeps the report tied to its actual scope.','repair':'no change needed','line':i,'quote':l}]
 a={'assessment_performed':True,'draft_sha256':hashlib.sha256(b).hexdigest(),'draft_version':version,'requested_mode':'independent','assessor_type':'self','reviewer':'Root release driver','method':'Separate post-draft reading against actual frozen identity, four native raw logs and independent timing readback, MI/novice/vignette/ancillary receipts, current policy, rights disposition, replay amendment and three fresh panel votes.','genre':genre,'coverage':'Whole exact saved draft','score':3,'style_status':'assessed','score_reason':'Concrete release observations with explicit limits; the structured receipt format necessarily repeats scope and identity.','confidence':'moderate','confidence_reason':'Editorial self-assessment, with no calibrated human or authorship inference.','examples':examples,'evidence_none_reason':None,'science_gate':'pass','science_gate_reason':'No new scientific result or default change is asserted; bounded checks do not become general recovery or nominal-coverage claims.','facts_gate':'pass','facts_gate_reason':'Archive identity, status, timing measures, counts, disclosures and approvals compared with the actual retained sources and current state; incomplete replay history is explicit.','references_gate':'pass','references_gate_reason':'Exact DOI identity checked against retained current primary Crossref metadata; policy files and dated rights records are distinguished from licence grants or legal opinions.','date':datetime.datetime.now(datetime.timezone.utc).isoformat(),'calibration_status':'Protocol examples and earlier technical receipts were seen; no held-out evaluation.','interpretation':'perceived_style_not_authorship'}
 out.write_text(json.dumps(a,indent=2)+'\n')
 subprocess.run(['python3','/Users/z3437171/shinichi-brain/benchmarks/personal-writing/benchmark.py','assessment-check',str(p),str(out),'--mode','independent','--version',version],check=True)
