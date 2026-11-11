"""核对三时期逐站降水与模板的字段。"""
from collections import defaultdict
from .io import date,number,positive

def summarize(v):
    counts=defaultdict(int);keys=set()
    for name in ('OH','MH','MF','template'):
        for r in v[name]:
            date(r['day']);key=(name,r['station'],r['day'])
            if key in keys:raise ValueError('同列站日重复')
            keys.add(key)
            if not isinstance(r['season'],str) or not r['season']:raise ValueError('季节身份缺失')
            if number(r['mm'])<0:raise ValueError('降水不得负')
            counts[name,r['station']]+=1
    if positive(v['wet_floor'])<=number(v['dry_threshold']):raise ValueError('湿日下限应高于干日阈值')
    return [{'period':p,'station':s,'records':n} for (p,s),n in sorted(counts.items())]
