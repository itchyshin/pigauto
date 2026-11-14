from pathlib import Path
import datetime,hashlib,json,time,urllib.request
r=Path(__file__).resolve().parent
identity=json.loads((r/'identity.json').read_text());archive=Path(identity['path']);digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert digest(archive)==identity['sha256']
url='https://api.crossref.org/works/10.1126/science.1157704';started=time.monotonic();request=urllib.request.Request(url,headers={'User-Agent':'pigauto-release-evidence/0.11 (read-only identifier verification)'})
with urllib.request.urlopen(request,timeout=20) as response:
 raw=response.read();status=response.status;final_url=response.url
(r/'crossref-raw.json').write_bytes(raw);m=json.loads(raw)['message']
expected_title='A Phylogenomic Study of Birds Reveals Their Evolutionary History'
assert status==200 and m['DOI']=='10.1126/science.1157704' and m['title']==[expected_title] and m['author'][0]['family']=='Hackett'
receipt={'url':url,'method':'GET','no_retry':True,'http_status':status,'final_url':final_url,'sha256':hashlib.sha256(raw).hexdigest(),'doi':m['DOI'],'title':m['title'],'journal':m['container-title'],'volume':m['volume'],'issue':m['issue'],'page':m['page'],'published':m['published'],'first_author':m['author'][0]['family'],'identity_matches':True,'elapsed_seconds':round(time.monotonic()-started,4),'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':digest(Path(__file__)),'archive_sha256_before':identity['sha256'],'archive_sha256_after':digest(archive),'doi_access':'Exact-source URL checker still returned403; metadata identity does not make that URL clean.'}
assert receipt['archive_sha256_after']==identity['sha256'];(r/'doi-metadata.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
