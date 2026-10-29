"""Small public-asset cache with exact bytes, durable failures, and bounded reads."""
from __future__ import annotations
import hashlib
import json
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[2]

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def utcnow():
    return datetime.now(timezone.utc).isoformat()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(path)

class PublicCache:
    def __init__(self, manifest):
        self.manifest = Path(manifest)
        self.assets = json.loads(self.manifest.read_text()) if self.manifest.exists() else []
        self._lock=threading.Lock()
        self._threads=threading.local()

    def _session(self):
        if not hasattr(self._threads,'session'):
            session=requests.Session()
            session.headers['User-Agent']='SpeaClimateResearch/0.1 (public-data feasibility)'
            retry=Retry(total=1,backoff_factor=1,status_forcelist=[429,500,502,503,504],respect_retry_after_header=True)
            session.mount('https://',HTTPAdapter(max_retries=retry))
            self._threads.session=session
        return self._threads.session

    def get_json(self, asset_id, url, params=None, source_id='S18', **metadata):
        path = ROOT / 'data/raw/gbif' / f'{asset_id}.json'
        with self._lock:
            prior = next((a for a in self.assets if a['asset_id'] == asset_id), None)
        if prior:
            if not path.exists() or sha256(path) != prior['sha256']:
                raise ValueError(f'Cached asset changed: {asset_id}')
            # Occurrence page size may adapt; offset/taxon and other query fields
            # must retain their identity. Use the returned row count to continue.
            for key,value in (params or {}).items():
                if key!='limit' and prior.get('request_parameters',{}).get(key)!=value:
                    raise ValueError(f'Cached query identity mismatch: {asset_id}/{key}')
            return json.loads(path.read_text())
        started=time.monotonic();received=0
        try:
            with self._session().get(url,params=params,timeout=(10,8),stream=True) as response:
                response.raise_for_status()
                chunks=[]
                for chunk in response.iter_content(8192):
                    received+=len(chunk)
                    if received>30_000_000:raise RuntimeError('Single public JSON response exceeded30MB; use smaller pages')
                    if time.monotonic()-started>45:raise TimeoutError('Public JSON response exceeded45seconds; use smaller pages')
                    chunks.append(chunk)
                payload=b''.join(chunks)
                result=json.loads(payload)
                response_url=response.url
                media_type=response.headers.get('Content-Type')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            asset=dict(asset_id=asset_id, source_id=source_id,
                citation='GBIF live API pilot; no download DOI',request_url=response_url,
                request_parameters=params or {},retrieved_at=utcnow(),
                release_version='live GBIF API v1',original_filename=path.name,
                media_type=media_type,sha256=sha256(path),byte_count=len(payload),
                license='Per-record licenses retained; metadata GBIF terms',
                geography='global queried; CA US MX audited consistently',time_coverage='all available',
                acquisition_status='fetched',path=str(path.relative_to(ROOT)),**metadata)
            with self._lock:
                self.assets.append(asset)
                write_json(self.manifest,self.assets)
            return result
        except (requests.RequestException,ValueError,TimeoutError,RuntimeError) as exc:
            failure={'asset_id':asset_id,'request_url':url,'request_parameters':params or {},
                     'retrieved_at':utcnow(),'status':'failed_attempt','exception_type':type(exc).__name__,
                     'elapsed_seconds':round(time.monotonic()-started,3),'partial_decoded_bytes_received':received,
                     'message':str(exc)[:500]}
            with self._lock:
                target=ROOT/'catalog/gbif_failures.jsonl'
                target.parent.mkdir(parents=True,exist_ok=True)
                with target.open('a') as handle:handle.write(json.dumps(failure)+'\n')
            raise
