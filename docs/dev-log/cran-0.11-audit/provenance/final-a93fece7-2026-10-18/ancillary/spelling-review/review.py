import csv,json,pathlib,hashlib,tarfile,collections,re
BASE=pathlib.Path('/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/ancillary')
PROOF=BASE; SRC=PROOF/'exact-source/pigauto'; OUT=BASE/'spelling-review'
ARCH=pathlib.Path('/private/tmp/pigauto-release-candidate/sha256-a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb/pigauto_0.11.0.tar.gz')
SHA='a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb'
assert hashlib.sha256(ARCH.read_bytes()).hexdigest()==SHA
rows=list(csv.DictReader((PROOF/'spelling-results-flat.csv').open())); assert rows and len({r['word'] for r in rows})==len(rows)
G={}
def group(k,s):
 for w in s.split():
  assert w not in G,w
  G[w]=k
group('accepted_British_variant','analysing artefact artefacts behaviour behavioural centre centred characterisation colour coloured factorisation favour generalised generalises honours initialisation initialisations initialised initialises labelled licence maximising minimise minimises minimising neighbourhood neighbours normalisation normalised optimised optimiser parallelising parameterised penalises personalised recognises regularisation regularisations regularised relabelling renormalises reorganisation reorganises serialised Serialises stabilising summarise Summarise summarises unnormalised visualisations Visualise')
group('mathematical_code_DOI_token','acc AE AO arg args basisOfRecord BCE bigl bigr BodyShapeI boldsymbol bw cbinds CDN cdot cffa CLR CMD compat concat coords cov covs CTmax CUDA cv DAE de decd DEP DepthRangeDeep dev df DGP diag doi enum env eqn expm et al ffn FFN FMI frac gbif GH GLS gnn GNN GNN\'s gzip hasCoordinate hasGeospatialIssues HTMLs ij ik iter iters JB js kk L̂ ldots LKJ lme LOC lon L̂ᵀR mathbb mathbf mathcal mathrm mathtt MCAR MCC MCSE md MLP MNAR MPS mult MVN nrow occ oo OO OOM OOM\'d optim otimes OU OVR PGLS phy phylo Phylo pmm PMM pre Pre preds propto qquad rds README RJ rmse rowMeans rsbl sd SDs SEs SHA sla SLA SLURM SMCFCS sysbio Syst syy th Troph un underbrace UQ vcv vec venv Vphy VU widehat xz zi ZI zRMSE Δ λ μ ρ σ Σ Σ̂')
group('valid_scientific_statistical_term','activations argmax autoencoder Autoencoder bioclim biogeographic biome BLUPs clade Clade clades cladogram compositional culmen denoising Denoising dimensionally discriminative eigen eigendecomposition eigensolver estimands exchangeability inferentially Insessorial interpretability latents learnable lossiness lossy microbiome mis mispaired misspecification morphometric morphometrics Nonresponse overfits patristic pectinate phylogenetically phylogenies Phylogenies phylogenomic polytomies polytomy poolable pretraining recomputation scaler softmax softplus stochasticity suboptimal subspaces terrestriality Topologies trophic Trophic ultrametric unconverged undercoverage underperform undersample undersampled')
group('valid_computing_ordinary_compound','allocator autoresearch Autoresearch benchmarked benchmarking comparator discoverability docstring downloader dropdown errored extractability fallbacks hardcoded lifecycle Lifecycle navbar passthrough precompute precomputes predeclared preflight preprocess Preprocess schemas selectable subprocess supersession transformative uncached validator worktree writeup')
group('proper_taxon_package_name','AdamW Aitchison AmNat AmphiBIO Ané AVONET BACE BACE\'s bennetti BES BIEN BIEN\'s Biometrika BirdTree BirdTree\'s Blomberg\'s Britton brms brmsfit Bruggeman Buerkner Buuren Casuarius Cohn Delhey drmTMB Ericson ESA Felsenstein Fick FishBase fishtree Freckleton Gammerman GBIF Gelman ggplot Ghahramani gllvmTMB glmmTMB GlobTherm Goolsby Grafen Graphormer Groothuis Hackett Hadfield Hartmann Heringa Hijmans IUCN Jetz jomo Kimball Krylov Lanczos MCMCglmm MEE mondrian Mooers MrBayes NAR Narval Neoaves NeurIPS Newick nlme Ornstein Oudshoorn Pagel Pagel\'s PAICE PanTHERIA Penone PhyloMaker phylopars PhyloPars phytools pigauto\'s pkgdown Rabosky Rao Reddy Reiter ReLU ResidualPhyloDAE ResNet rfishbase\'s rgbif roxygen Roxygen Rphylopars RSpectra Shafer Shinichi\'s Sociobiology SoilGrids TabPFN TabPFN\'s testthat Tobias Totoro Uhlenbeck Vaswani Vehtari Villemereuil Vovk WorldClim Ying Zhou')
unknown=[r['word'] for r in rows if r['word'] not in G]
assert not unknown,unknown
contexts=[]; result=[]; used={}; issues=[]; cache={}
with tarfile.open(ARCH,'r:gz') as tar:
 for row in rows:
  word=row['word']; scopes=set(); matches=[]
  for ref in row['found'].split(';'):
   name,ns=ref.strip().rsplit(':',1)
   path=SRC/name
   if name=='description':path=SRC/'DESCRIPTION'
   elif not path.exists():path=SRC/('man' if name.endswith('.Rd') else 'vignettes')/name
   assert path.is_file(),path
   rel=str(path.relative_to(SRC)); data=path.read_bytes(); member=cache.setdefault(rel, tar.extractfile('pigauto/'+rel).read()) if rel not in cache else cache[rel]
   assert data==member,rel
   used[rel]=hashlib.sha256(data).hexdigest(); lines=data.decode().splitlines()
   scope='historical_NEWS' if name=='NEWS.md' else 'current_reader_documentation'
   scopes.add(scope)
   for n in map(int,ns.split(',')):
    candidates=[i for i,line in enumerate(lines,1) if any((m.start()==0 or not line[m.start()-1].isalpha()) and (m.end()==len(line) or not line[m.end()].isalpha()) for m in re.finditer(re.escape(word),line))]
    assert candidates,(word,name)
    # spelling locations may identify the enclosing paragraph; retain the closest actual token line.
    actual=min(candidates,key=lambda i:abs(i-n))
    context=' | '.join(lines[max(0,actual-2):min(len(lines),actual+1)])
    contexts.append({'word':word,'reported_reference':name+':'+str(n),'source_path':rel,'actual_line':actual,'scope':scope,'context':context})
    matches.append(f'{rel}:{actual}')
  cat=G[word]; note={
  'accepted_British_variant':'Accepted British spelling; no spelling correction indicated.',
  'mathematical_code_DOI_token':'Identifier, abbreviation, notation, reference particle, URL/DOI component or token fragment; not an ordinary-word typo.',
  'valid_scientific_statistical_term':'Valid technical term or productive technical compound in this context; no spelling correction indicated.',
  'valid_computing_ordinary_compound':'Valid computing term, productive compound, or informal historical wording; no spelling correction indicated.',
  'proper_taxon_package_name':'Name, named method, taxon, package, journal or exact API value; proper-token classification does not verify bibliographic identity.'}[cat]
  if word in ['cffa','decd']:note='Speller fragment of hexadecimal commit identifier ('+('dc8cffa' if word=='cffa' else '573decd')+'); no word correction indicated.'
  if word in ['bigl','bigr']:note='LaTeX delimiter-sizing command; no word correction indicated.'
  if word=='mis':note='Productive prefix in mis-extrapolates and mis-pool; no standalone prose typo.'
  if word=='un':note='Productive prefix in un-log/un-restored; wording is technical shorthand, not a misspelling.'
  if word=='arg':note='Historical NEWS shorthand for argument; expansion is optional editorial style, not an ordinary-word typo.'
  if word=='passthrough':note='Historical NEWS computing compound for passing an input unchanged; no word correction indicated.'
  if word=='mondrian':note='Exact quoted conformal_method API value; do not capitalize inside code.'
  if word=='writeup':note='Informal historical NEWS compound (often styled write-up); optional editorial hyphenation, not an ordinary-word typo.'
  result.append({'word':word,'category':cat,'adjudication':'no_spelling_change_indicated','scope':';'.join(sorted(scopes)),'reported_references':row['found'],'resolved_contexts':';'.join(matches),'note':note})
for name,data in [('token-adjudication.csv',result),('source-contexts.csv',contexts)]:
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
summary={'archive_sha256':SHA,'archive_bytes':ARCH.stat().st_size,'input_sha256':hashlib.sha256((PROOF/'spelling-results-flat.csv').read_bytes()).hexdigest(),'flag_rows':len(rows),'unique_tokens':len(set(r['word'] for r in rows)),'resolved_context_rows':len(contexts),'source_files_bound_to_archive':used,'categories':dict(collections.Counter(r['category'] for r in result)),'scopes':dict(collections.Counter(r['scope'] for r in result)),'ordinary_word_typos_identified':[],'uncertain_flags':[],'limitations':['Triage covers only the actual fresh flagged tokens, not all unflagged prose.','Names and bibliographic fragments were assessed as spelling tokens, not independently reference-verified.','Spelling checker line references can mark a paragraph rather than the exact token; source-contexts.csv retains original and nearest actual lines.','No dictionary, package source, archive, workflow or warning audit was edited.']}
(OUT/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='source_files_bound_to_archive'},indent=2))
