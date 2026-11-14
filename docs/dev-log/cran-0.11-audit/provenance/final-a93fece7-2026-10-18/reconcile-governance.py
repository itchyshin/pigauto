"""Advance current records to the proven successor; preserve all prior records."""
from pathlib import Path
import datetime,hashlib,json,shutil
S=Path(__file__).resolve().parent;R=Path('/private/tmp/pigauto-pr236-validator-20261017');A=R/'docs/dev-log/cran-0.11-audit';O=A/'provenance/final-a93fece7-2026-10-18';P=str(O.relative_to(R))+'/'
F=json.loads((S/'FINAL.json').read_text());old=json.loads((A/'release-ledger.json').read_text())
assert F['sha256']=='a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb'
assert json.loads((O/'actual-warning-audit/verdict.json').read_text())['verdict']=='PASS'
assert json.loads((O/'actual-warning-audit/verdict.json').read_text())['passed_cells']==4
assert not (A/'historical-release-ledger-pre-a93fece7-2026-10-18.json').exists()
shutil.copyfile(A/'release-ledger.json',A/'historical-release-ledger-pre-a93fece7-2026-10-18.json')
shutil.copyfile(R/'cran-comments.md',S/'cran-comments-before-successor.md')
shutil.copyfile(S/'cran-comments-draft.md',R/'cran-comments.md')
def evidence(p):
 assert (O/p).is_file(),p
 return {'kind':'file','value':P+p}
persistent=Path('/Users/z3437171/.codex/artifacts/pigauto-cran-0.11.0')/('sha256-'+F['sha256'])/'pigauto_0.11.0.tar.gz'
assert hashlib.sha256(persistent.read_bytes()).hexdigest()==F['sha256']
x={
'package':'pigauto','version':'0.11.0','release_type':'update','status_claim':'platform-clean','next_unproven_rung':'submission-ready',
'audit_verdict':'NOT_READY_PENDING_FRESH_FINAL_PANEL_AND_EXECUTABLE_LEDGER',
'artifact':{k:F[k] for k in ['sha256','size_bytes','source_commit','clean_worktree','inventory_sha256','archive_entries','forbidden_entries']},
'profile':old['profile'],'evidence':{},'panel':{n:'NOT_READY' for n in ['grace','rose','pat']},'panel_evidence':{},
'submission':old['submission'],'evidence_pr':old['evidence_pr'],
'source_candidate':{'source_commit':F['source_commit'],'source_pr':238,'source_pr_merged':True,'approved_pr_head':'816987fd160c1d2c5c58b020b07357ff17c077e5','main_source_ci':39090594763,'source_ci_state':'PASS_ALL3_MERGED_HEAD','pages_run':39090594794,'deployed_site_gate':'PASS'},
'approvals':{'source_merge_and_automatic_deployment':'Explicit Merge238 approval exercised and actual deployment verified at b2e6054','public_builder_upload':'Explicitly authorized; existing validated R-hub account used, no new token request; public original gzip fullhash/size match','cran_upload':'Conditional approval after all gates pass; not submitted','maintainer_confirmation':'Included in conditional submission authority; not performed','evidence_pr236_merge':'Not authorized; remain open draft unmerged','tag_and_announcement':'Not authorized; not performed'},
'historical_ledger':'docs/dev-log/cran-0.11-audit/historical-release-ledger-pre-a93fece7-2026-10-18.json',
'predecessor':{'sha256':old['artifact']['sha256'],'state':'withdrawn and immutable; old panel votes not promoted'},
'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
'claim_scope':'Exact original gzip native checks on four platforms, actual warning-content review, installed optionalmodel/novice and ancillary diagnostics; final independent release panel pending.',
'replay_provenance_limit':'Two successful supporting captures retained with separate coherent timing pairs; prior mixed derived timing corrected by amendment. Identical warning-object bytes match both, so unique capture attribution is unresolved. Native checks were not rerun.',
'final_proof_summary':P+'README.md',
}
x['artifact'].update(path=str(persistent),original_checked_path=F['path'],inventory_evidence=evidence('inventory.txt') if (O/'inventory.txt').exists() else evidence('freeze-review/actual-verdict.json'),current_candidate=True,state='Frozen successor; four native checks StatusOK; final panel pending')
for k,p in {'product_contract':'README.md','rights_and_consent':'README.md','current_cran_policy':'policy-refresh/receipt.json','rendered_site':'site/verdict.json','local_as_cran':'native/macos-arm64/Rcheck/00check.log','incoming_feasibility':'native-platform-summary.json','installed_package':'installed-novice/verdict.json','timing':'native-platform-summary.json','platform_matrix':'native/direct-verification.log','external_logs':'native/server-digest-binding.json','external_service_behaviour':'ancillary/verdict.json','size_and_vignette_budget':'vignette-timing/verdict.json'}.items():
 # README is generated next; all other concrete evidence is required now.
 if p=='README.md':x['evidence'][k]={'kind':'file','value':P+p}
 else:x['evidence'][k]=evidence(p)
x['evidence']['reverse_dependencies']={'kind':'rationale','value':'NOT_APPLICABLE: primary CRAN PACKAGES.gz read on 2026-10-18 measured zero reverse dependencies across Depends, Imports, LinkingTo, Suggests and Enhances; exact source receipt '+P+'cran-publication-and-revdeps.json'}
x['evidence']['compiled_diagnostics']={'kind':'rationale','value':'NOT_APPLICABLE: frozen inventory and source-size binding show no package-owned src or compiled binaries; dependency runtimes were tested separately.'}
x['timing_decision']={'native_status':'OK on all4, noNOTE','native_timing_evidence':P+'native-platform-summary.json','vignette_rebuilds':'15secondsMac elapsed,20secondsWindowsrelease,17secondsWindowsdevel,19secondsLinux elapsed','separate_replay_timing':'Two pairs disclosed in actual-warning-audit/timing-provenance-amendment; excluded from native check duration','independent_final_disposition':'PENDING_FRESH_PANEL'}
(A/'release-ledger.json').write_text(json.dumps(x,indent=2)+'\n')
shutil.copyfile(S/'persistent-readback.json',O/'persistent-readback.json')
shutil.copyfile(S/'publish-successor-evidence.py',O/'publish-successor-evidence.py')
shutil.copyfile(S/'reconcile-governance.py',O/'reconcile-governance.py')
print('SUCCESSOR_GOVERNANCE_RECONCILED_PLATFORM_CLEAN_PANEL_PENDING')
