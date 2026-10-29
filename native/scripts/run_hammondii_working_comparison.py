"""Bounded empirical conditional-density comparison; declared assumptions, no tests of significance."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key] = '1'
import argparse
import json
import resource
import sys
import time
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.special import logsumexp
from spea.provenance import ROOT,sha256,utcnow,write_json
from spea.fixed_basis_density import fit_conditional_basis
from spea.smooth_density_basis import SmoothDensityBasis
from spea.genus_model import geographic_blocks,stratified_quadrature


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--plan',default='config/hammondii_working_comparison.json')
    args=ap.parse_args()
    start,cpu=time.monotonic(),time.process_time()
    out=ROOT/'runs'/args.run_id;out.mkdir(exist_ok=False)
    for name in ['splits','fits','full_predictions']:(out/name).mkdir()
    cfg=json.loads((ROOT/args.plan).read_text());write_json(out/'effective_config.json',cfg)
    inputs,checks,metrics,scores,audit,coverage=[],[],[],[],[],[]
    hashes={}
    def check(name,okay):
        checks.append(dict(name=name,passed=bool(okay)))
        if not okay:raise ValueError(name)
    def asset(run,name):
        directory=ROOT/'runs'/run;manifest=json.loads((directory/'run.json').read_text());p=directory/name
        check('hash:'+str(p.relative_to(ROOT)),sha256(p)==manifest['output_hashes'][name])
        for q in [p,directory/'run.json']:hashes[str(q.relative_to(ROOT))]=sha256(q)
        return p
    def status(value,**kw):
        write_json(out/'run.json',dict(status=value,created_utc=utcnow(),command_argv=sys.argv,config=cfg,
            new_fits=len(audit),empirical_fits=len(audit),formal_significance_tests=0,
            cpu_seconds=time.process_time()-cpu,elapsed_seconds=time.monotonic()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**kw))
    def flush():
        for name,rows in [('fit_audit',audit),('metrics',metrics),('cell_scores',scores),('coverage',coverage)]:
            if rows:pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    status('preparing')
    grid=dict(np.load(asset(cfg['geography_run'],'full_geography.npz')))
    available=dict(np.load(asset(cfg['support_run'],'available_geography.npz')))
    cells=pd.read_csv(asset(cfg['support_run'],'cell_geographic_assignments.csv'))
    cohort=np.load(asset(cfg['source_audit_run'],'cohort_native_ids.npz'))
    reg=json.loads(asset(cfg['preparation_run'],'feature_registry.json').read_text())
    # The existing audit separates linked source/event conflicts from calibration edges.
    role=pd.read_csv(asset('2026-09-27-paired-hammondii-position-audit-v1','cell_fold_diagnostics.csv.gz'))
    native,g,area,block=[grid[k] for k in ['native','group','area_m2','block50']]
    cp=np.searchsorted(native,cells.native_flat.to_numpy())
    check('all current presence cells in full geometry',np.array_equal(native[cp],cells.native_flat))
    check('unchanged groups',np.array_equal(g[cp],cells.provisional_group.map({'North':1,'South':-1})))
    aptr=np.searchsorted(available['native'],native)
    check('same full availability',np.array_equal(available['native'][aptr],native))
    distance=available['signed_distance_m'][aptr]
    fields={f['key']:f for f in reg['fields']};climate=[]
    for key in cfg['features']:
        p=ROOT/fields[key]['physical_path'];asset(p.parent.name,p.name)
        v=np.load(p,mmap_mode='r').reshape(-1)[native].astype(float)
        if key=='bio12':
            check('nonnegative rainfall',np.all(v>=0));v=np.log1p(v)
        climate.append(v)
    x=np.column_stack(climate);check('all climate finite',np.isfinite(x).all())
    no35=np.isin(native,np.asarray(cohort['without_S35']))
    check('primary cohort equals current cells',np.array_equal(np.sort(cohort['primary']),np.sort(native[cp])))
    check('source subset',np.isin(cohort['without_S35'],native[cp]).all())
    np.savez_compressed(out/'analysis_grid.npz',native=native,group=g,area_m2=area,block50=block,
        distance_m=distance,climate=x,presence_index=cp,without_S35=no35,
        longitude=grid['longitude'],latitude=grid['latitude'])
    block25=geographic_blocks(native,reg['window_grid']['shape'],reg['window_grid']['transform'],25000)
    source=[Path(__file__),ROOT/args.plan,ROOT/'src/spea/smooth_density_basis.py',ROOT/'src/spea/fixed_basis_density.py',
            ROOT/'src/spea/genus_model.py',ROOT/'uv.lock',ROOT/'pyproject.toml']
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in source:z.write(p,str(p.relative_to(ROOT)));hashes[str(p.relative_to(ROOT))]=sha256(p)
    geod=Geod(ellps='WGS84');split_cache={};basis_cache={}
    def prepare(case,fold):
        key=f'{case["name"]}_{fold}'
        if key in split_cache:return split_cache[key]
        full=fold=='full';held=np.array([],dtype=int) if full else cp[grid[f'held{fold}'][cp]]
        domain=np.arange(len(native)) if full else np.flatnonzero(grid[f'train{fold}'])
        train_p=cp if full else cp[grid[f'train{fold}'][cp]]
        radii=cells[f'radius_unknown{case["unknown_radius_km"]}km_m'].to_numpy()
        radius_by_grid=np.zeros(len(native));radius_by_grid[cp]=radii
        # Exact pair guard, including both occurrence supports; no finite-M containment filter.
        if not full:
            a,b=np.meshgrid(train_p,held,indexing='ij')
            sep=np.abs(geod.inv(grid['longitude'][a],grid['latitude'][a],grid['longitude'][b],grid['latitude'][b])[2])
            bound=sep-radius_by_grid[a]-radius_by_grid[b]
            eligible=bound.min(axis=1)>=30000
            train_p=train_p[eligible]
            conflicts=role[(role.fold==fold)&(role.scenario==f'{case["unknown_radius_km"]}km')]
            train_p=train_p[~np.isin(native[train_p],conflicts.loc[conflicts.source_guard_conflict,'native_flat'])]
            check(key+':disjoint train/held blocks',not np.isin(block[train_p],block[held]).any())
        forced=train_p.copy()
        if case['source']=='without_S35':train_p=train_p[no35[train_p]]
        qids,qw,pi=[],[],[]
        for gg in [-1,1]:
            dd=domain[g[domain]==gg];ff=forced[g[forced]==gg]
            q,w,p,_,_=stratified_quadrature(native[dd],native[ff],block25[dd],area[dd],
                cfg['quadrature']['nonpresence_per_25km_block'],cfg['quadrature']['seed']+(0 if full else int(fold))+gg+1)
            qids.append(np.searchsorted(native,q));qw.append(w);pi.append(p)
        qids,qw,pi=map(np.concatenate,[qids,qw,pi]);order=np.argsort(qids);qids,qw,pi=qids[order],qw[order],pi[order]
        check(key+':presences forced into quadrature',np.isin(train_p,qids).all())
        target=np.arange(len(native)) if full else np.flatnonzero(grid[f'held{fold}'])
        mass=np.zeros(len(qids))
        for gg in [-1,1]:
            pp=train_p[g[train_p]==gg]
            check(key+f':train minimum {gg}',len(pp)>=cfg['minimum_train_cells_each_group'])
            mass[np.searchsorted(qids,pp)]=.5/len(pp)
            coverage.append(dict(case=case['name'],fold=str(fold),group=gg,train_cells=len(pp),
                train_blocks=len(np.unique(block[pp])),held_cells=int((g[held]==gg).sum()),
                held_blocks=len(np.unique(block[held[g[held]==gg]])),
                held_without_S35=int(((g[held]==gg)&no35[held]).sum()),quadrature_cells=int((g[qids]==gg).sum())))
        result=dict(q=qids,w=qw,pi=pi,mass=mass,p=train_p,held=held,target=target,domain=domain)
        np.savez_compressed(out/'splits'/(key+'.npz'),**result)
        split_cache[key]=result
        return result
    # Freeze geometric/source inputs before any fit.
    tasks=[]
    for case in cfg['cases']:
        for fold in cfg['folds']:
            prepare(case,fold)
            for geometry in ['spline']+(['linear'] if case['name']=='all_1km' else []):
                for ridge in cfg['ridges']:
                    for structure in cfg['structures']:tasks.append((case,fold,geometry,ridge,structure))
    for case in cfg['cases'][:2]:
        prepare(case,'full')
        for ridge in cfg['ridges']:
            for structure in cfg['structures']:tasks.append((case,'full','spline',ridge,structure))
    check('frozen fit count',len(tasks)==cfg['budget']['maximum_fits'])
    write_json(out/'prefit_receipt.json',dict(input_hashes=hashes,planned_fits=len(tasks),new_fits=0,
        frozen_outputs={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file() and p.name!='run.json'}))
    def predict(basis,coeff,cols,indices):
        result=np.empty(len(indices))
        for k in range(0,len(indices),cfg['chunk_size']):
            ii=indices[k:k+cfg['chunk_size']]
            result[k:k+len(ii)]=basis.transform(x[ii],distance[ii],g[ii])[:,cols]@coeff
        return result
    try:
        for case,fold,geometry,ridge,structure in tasks:
            if time.process_time()-cpu>cfg['budget']['cpu_minutes']*60:raise RuntimeError('CPU budget exceeded')
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>cfg['budget']['peak_gib']*2**30:raise RuntimeError('Memory budget exceeded')
            key=f'{case["name"]}_{fold}';split=split_cache[key];q,w,mass=split['q'],split['w'],split['mass']
            bkey=key+'_'+geometry
            if bkey not in basis_cache:
                basis=SmoothDensityBasis().fit(x[q],distance[q],g[q],w,geography=geometry)
                design=basis.transform(x[q],distance[q],g[q]);basis_cache[bkey]=(basis,design)
                write_json(out/'splits'/(bkey+'_basis.json'),basis.metadata())
            basis,design=basis_cache[bkey];cols=basis.columns[structure]
            before=time.process_time()
            fit=fit_conditional_basis(design[:,cols],g[q],w,mass,ridge=ridge,tol=cfg['optimization_tolerance'])
            stem=f'{key}_{geometry}_{structure}_ridge{ridge:g}'
            np.savez_compressed(out/'fits'/(stem+'.npz'),**fit)
            target,held=split['target'],split['held']
            eta=predict(basis,fit['coefficients'],cols,target)
            norm={gg:float(logsumexp(np.log(area[target[g[target]==gg]])+eta[g[target]==gg])-np.log(area[target[g[target]==gg]].sum())) for gg in [-1,1]}
            normalized=eta-np.array([norm[int(gg)] for gg in g[target]])
            for gg in [-1,1]:
                use=g[target]==gg
                check(stem+f':area normalization {gg}',np.isclose(np.sum(area[target[use]]/area[target[use]].sum()*np.exp(normalized[use])),1,atol=1e-10))
            if fold=='full':
                np.savez_compressed(out/'full_predictions'/(stem+'.npz'),log_relative_density=normalized,log_intensity=eta)
            else:
                ti=np.searchsorted(target,held);check(stem+':held containment',np.array_equal(target[ti],held))
                y=normalized[ti]
                for i,ii in enumerate(held):
                    scores.append(dict(model=stem,case=case['name'],fold=fold,geography=geometry,structure=structure,ridge=ridge,
                        native_flat=int(native[ii]),group=int(g[ii]),block50=int(block[ii]),without_S35=bool(no35[ii]),score=float(y[i])))
                for target_cohort,keep in [('all',np.ones(len(held),bool)),('without_S35',no35[held])]:
                    for gg in [-1,1]:
                        use=keep&(g[held]==gg);check(stem+f':held minimum {target_cohort}/{gg}',int(use.sum())>=cfg['minimum_test_cells_each_group'])
                        yy=y[use];bb=block[held[use]]
                        block_mean=pd.DataFrame(dict(score=yy,block=bb)).groupby('block').score.mean().mean()
                        metrics.append(dict(model=stem,case=case['name'],fold=fold,geography=geometry,structure=structure,ridge=ridge,
                            target_cohort=target_cohort,group=gg,cells=int(use.sum()),blocks=len(np.unique(bb)),
                            cell_score=float(yy.mean()),block_score=float(block_mean)))
            audit.append(dict(model=stem,case=case['name'],fold=str(fold),geography=geometry,structure=structure,ridge=ridge,
                columns=len(cols),objective=fit['objective'],gradient_inf=fit['gradient_inf'],iterations=fit['iterations'],
                fit_and_prediction_cpu_seconds=time.process_time()-before,north_normalizer=norm[1],south_normalizer=norm[-1],
                sha256=sha256(out/'fits'/(stem+'.npz'))))
            if len(audit)%8==0:
                flush();status('running');print(json.dumps(dict(fits=len(audit),cpu_seconds=time.process_time()-cpu)),flush=True)
        flush();check('all planned fits completed',len(audit)==len(tasks))
        for path,h in hashes.items():check('unchanged:'+path,sha256(ROOT/path)==h)
        write_json(out/'checks.json',dict(passed=len(checks),failed=0,checks=checks))
        status('completed',input_hashes=hashes,output_hashes={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file() and p.name!='run.json'})
    except Exception as error:
        flush();status('failed',error=repr(error),input_hashes=hashes);raise


if __name__=='__main__':main()
