"""Offline, interactive StrainHub-style country networks for HTML reports.

Only the presentation changes: displayed metrics are those of the representative
history, and possible alternative edges never replace that history. vis-network
9.1.13 is vendored from the official npm package under its upstream licences.
"""

from __future__ import annotations

import hashlib
import html
import json
import math
import re
from importlib.resources import files
from uuid import uuid4

from .colours import country_palette

_METRICS = (
    ("in_degree", "Indegree centrality"),
    ("out_degree", "Outdegree centrality"),
    ("degree", "Degree centrality"),
    ("betweenness", "Betweenness centrality"),
    ("closeness", "Closeness centrality"),
    ("source_hub_ratio", "Source hub ratio"),
)


def _number(value, default=0):
    """Keep non-finite numbers out of JSON and the layout engine."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else None


def _payload(row, *, focus_inputs=False, include_possible=False):
    metrics = {str(item["country"]): item for item in row.get("network_metrics", [])}
    directed = row.get("directed_edges", row.get("edges", []))
    originals = row.get("nodes", [])
    if not originals:
        # Older report fixtures contain endpoints but no explicit node table.
        countries = sorted({str(edge[key]) for edge in directed for key in ("source", "target")})
        originals = [
            dict(
                country=country, count=0, is_input_country=country in row.get("input_countries", [])
            )
            for country in countries
        ]
    palette = country_palette(str(node["country"]) for node in originals)
    palette.update(
        {
            str(key): value
            for key, value in row.get("country_colors", {}).items()
            if isinstance(value, str) and re.fullmatch(r"#[0-9a-fA-F]{6}", value)
        }
    )
    nodes = []
    for node in originals:
        country = str(node["country"])
        values = metrics.get(country, {})
        measures = {key: _number(values.get(key)) for key, _ in _METRICS}
        if measures["degree"] is None:
            measures["degree"] = (measures["in_degree"] or 0) + (measures["out_degree"] or 0)
        nodes.append(
            dict(
                id=country,
                country=country,
                color=palette[country],
                input=bool(node.get("is_input_country")),
                count=_number(node.get("count", 0)),
                metrics=measures,
            )
        )
    countries = {node["id"] for node in nodes}
    inputs = {node["id"] for node in nodes if node["input"]}
    edges = []
    for item in directed:
        source, target = str(item["source"]), str(item["target"])
        if source not in countries or target not in countries or source == target:
            continue
        count = _number(item.get("representative_count", 0)) or 0
        minimum = _number(item.get("min_changes", count)) or 0
        maximum = _number(item.get("max_changes", count)) or 0
        if max(count, minimum, maximum) <= 0:
            continue
        if focus_inputs and not ({source, target} & inputs):
            continue
        edges.append(
            dict(
                id=f"edge-{len(edges)}",
                source=source,
                target=target,
                count=count,
                minimum=minimum,
                maximum=maximum,
                ambiguous=bool(item.get("count_ambiguous", minimum != maximum)),
            )
        )
    if focus_inputs:
        retained = inputs | {edge[key] for edge in edges for key in ("source", "target")}
        nodes = [node for node in nodes if node["id"] in retained]
    return dict(nodes=nodes, edges=edges, includePossible=bool(include_possible))


def interactive_network_html(row, identifier, *, focus_inputs=False, include_possible=False):
    """Return a replaceable widget host; include :func:`widget_assets` once per report."""
    # A random suffix prevents collisions when one cohort occurs in multiple scopes.
    token = hashlib.sha256(str(identifier).encode()).hexdigest()[:12] + "-" + uuid4().hex[:12]
    payload = json.dumps(
        _payload(row, focus_inputs=focus_inputs, include_possible=include_possible),
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
    )
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    options = "".join(
        f'<option value="{key}">{html.escape(label)}</option>' for key, label in _METRICS
    )
    return f"""<section class="cc-network-widget" data-cc-network data-widget="country-network">
  <div class="cc-network-controls">
    <label for="cc-metric-{token}">Node size <select id="cc-metric-{token}" data-network-metric>{options}</select></label>
    <label for="cc-country-{token}">Highlight country <select id="cc-country-{token}" data-network-country><option value="">All countries</option></select></label>
    <label><input type="checkbox" data-network-possible{" checked" if include_possible else ""}> Include alternative possible links</label>
    <div class="cc-network-buttons" aria-label="Network navigation">
      <button type="button" data-network-action="zoom-in" aria-label="Zoom in">+</button>
      <button type="button" data-network-action="zoom-out" aria-label="Zoom out">−</button>
      <button type="button" data-network-action="fit">Fit network</button>
      <button type="button" data-network-action="clear">Clear selection</button>
    </div>
  </div>
  <p class="cc-network-metric-note" data-network-metric-note>Node size: indegree centrality.</p>
  <div class="cc-network-canvas" data-network-canvas role="img" aria-label="Interactive directed country network. Use the country selector and navigation buttons to explore connections."></div>
  <p class="cc-network-status" data-network-status role="status" aria-live="polite"></p>
  <div class="cc-network-legend"><span>Arrow width: representative change count</span><span>Dashed: count varies across optimal histories</span><span>Dark outline: input country</span></div>
  <details class="cc-network-country-key"><summary>Country colour key</summary><div class="cc-network-countries" data-network-legend aria-label="Country colours"></div></details>
  <p class="cc-network-caveat">Arrows show country-state changes in one optimal history, not demonstrated transmission. Alternative links need not occur together. Metrics remain those of the representative history.</p>
  <script type="application/json" data-network-data>{payload}</script>
</section>"""


_CSS = r"""
.cc-network-widget{border:1px solid var(--line,#d8dee6);border-radius:0;padding:14px;background:#fff;color:var(--ink,#263142);margin:12px 0}
.cc-network-controls{display:flex;align-items:center;flex-wrap:wrap;gap:12px;font-size:13px}
.cc-network-controls label{display:flex;align-items:center;gap:6px}.cc-network-controls select{max-width:230px}
.cc-network-controls select,.cc-network-buttons button{font:inherit;border:1px solid var(--line,#b8c2d1);border-radius:0;background:#fff;color:var(--ink,#263142);padding:6px 8px}
.cc-network-controls select:focus-visible,.cc-network-buttons button:focus-visible{outline:3px solid #ffb800;outline-offset:2px}
.cc-network-buttons{display:flex;gap:5px;flex-wrap:wrap}.cc-network-buttons button{cursor:pointer}
.cc-network-canvas{height:560px;min-height:340px;width:100%;border:1px solid var(--line,#e6eaf0);border-radius:0;background:#fafbfd}
.cc-network-widget .vis-tooltip{position:absolute;visibility:hidden;z-index:10;white-space:pre-line;padding:8px 10px;background:#fff;color:var(--ink,#263142);border:1px solid var(--line,#b8c2d1);font:12px/1.5 Arial,sans-serif;pointer-events:none;max-width:320px}
.cc-network-metric-note,.cc-network-status{margin:10px 0;font-size:13px}.cc-network-metric-note{font-weight:600}
.cc-network-legend,.cc-network-countries{display:flex;flex-wrap:wrap;gap:8px 18px;font-size:12px;margin:10px 0}
.cc-network-countries span{display:inline-flex;align-items:center;gap:5px}.cc-network-swatch{width:10px;height:10px;border-radius:50%;display:inline-block}
.cc-network-country-key{font-size:12px;margin-top:8px}.cc-network-country-key summary{cursor:pointer}
.cc-network-caveat{font-size:12px;line-height:1.5;color:#596579;margin:10px 0 0}
@media(max-width:600px){.cc-network-canvas{height:440px}.cc-network-controls{align-items:flex-start}.cc-network-controls label{flex-wrap:wrap}}
"""

_INITIALIZER = r"""
(function(){
  if(window.ChronoCladeNetworks) return;
  const instances=new Map();
  const names={in_degree:'Indegree centrality',out_degree:'Outdegree centrality',degree:'Degree centrality',betweenness:'Betweenness centrality',closeness:'Closeness centrality',source_hub_ratio:'Source hub ratio'};
  const visible=host=>host.isConnected && !host.closest('[hidden]') && host.getBoundingClientRect().width>0;
  function destroy(container){
    for(const [host,instance] of instances){
      if(!container || host===container || container.contains(host) || !host.isConnected){
        instance.abort.abort(); instance.network.destroy(); instances.delete(host);
      }
    }
  }
  function mount(container=document){
    for(const [host,instance] of instances){
      if(!visible(host)){instance.abort.abort();instance.network.destroy();instances.delete(host);}
    }
    const hosts=[...(container.matches && container.matches('[data-cc-network]')?[container]:[]),...container.querySelectorAll('[data-cc-network]')];
    for(const host of hosts){
      if(!visible(host) || instances.has(host)) continue;
      const status=host.querySelector('[data-network-status]');
      if(!window.vis || !window.vis.Network){status.textContent='The network viewer is unavailable.';continue;}
      const data=JSON.parse(host.querySelector('[data-network-data]').textContent);
      const metric=host.querySelector('[data-network-metric]'),country=host.querySelector('[data-network-country]');
      const possible=host.querySelector('[data-network-possible]'),note=host.querySelector('[data-network-metric-note]');
      const canvas=host.querySelector('[data-network-canvas]'),legend=host.querySelector('[data-network-legend]');
      country.replaceChildren(new Option('All countries',''));legend.replaceChildren();
      for(const node of data.nodes){
        country.add(new Option(node.country,node.id));
        const entry=document.createElement('span'),swatch=document.createElement('i');
        swatch.className='cc-network-swatch';swatch.style.backgroundColor=node.color;
        entry.append(swatch,document.createTextNode(node.country));legend.append(entry);
      }
      const nodes=new vis.DataSet(),edges=new vis.DataSet();
      const network=new vis.Network(canvas,{nodes,edges},{
        layout:{randomSeed:1729,improvedLayout:true},
        physics:{solver:'repulsion',repulsion:{nodeDistance:180,centralGravity:0.12,springLength:180,springConstant:0.03,damping:0.12},stabilization:{iterations:250,fit:true}},
        interaction:{dragNodes:true,dragView:true,zoomView:true,hover:true,keyboard:{enabled:false},navigationButtons:false},
        nodes:{shape:'dot',font:{size:20,color:'#263142',face:'Arial'},scaling:{min:12,max:38,label:{drawThreshold:0}},borderWidth:1},
        edges:{arrows:{to:{enabled:true,scaleFactor:0.75}},smooth:{enabled:true,type:'dynamic'},color:{color:'#78869a',highlight:'#35445d',opacity:0.7},scaling:{min:1,max:8}}
      });
      const abort=new AbortController();instances.set(host,{network,abort});
      const listen=(element,event,callback)=>element.addEventListener(event,callback,{signal:abort.signal});
      function fitNetwork(){network.stopSimulation();network.fit({animation:false});network.moveTo({scale:network.getScale()*0.85,animation:false});}
      let activeEdges=[];
      function tooltip(text){const element=document.createElement('div');element.textContent=text;return element;}
      const format=value=>value===null?'undefined':Number(value.toPrecision(3)).toString();
      function update(){
        const selected=country.value,key=metric.value;
        activeEdges=data.edges.filter(edge=>edge.count>0 || (possible.checked && Math.max(edge.minimum,edge.maximum)>0));
        const neighbours=new Set(selected?[selected]:[]);
        for(const edge of activeEdges){if(edge.source===selected) neighbours.add(edge.target);if(edge.target===selected) neighbours.add(edge.source);}
        const numeric=data.nodes.map(node=>node.metrics[key]).filter(value=>value!==null && Number.isFinite(value));
        const maximum=Math.max(0,...numeric);
        nodes.update(data.nodes.map(node=>{
          const value=node.metrics[key],dim=selected && !neighbours.has(node.id);
          const background=dim?'#e2e6ed':node.color,border=node.input?'#222222':'#ffffff';
          return {id:node.id,label:node.country,size:12+(maximum>0 && value!==null?26*Math.sqrt(Math.max(0,value)/maximum):0),
            color:{background,border,highlight:{background,border},hover:{background,border}},borderWidth:node.input?3:1,
            font:{color:dim?'#a4acb9':'#263142'},title:tooltip(node.country+'\n'+names[key]+': '+format(value)+'\nSamples: '+node.count)};
        }));
        edges.clear();edges.add(activeEdges.map(edge=>({
          id:edge.id,from:edge.source,to:edge.target,value:Math.max(1,edge.count),dashes:edge.ambiguous || edge.count===0,
          color:{color:selected && !(edge.source===selected || edge.target===selected)?'#e2e6ed':'#78869a',highlight:'#35445d',opacity:0.7},
          title:tooltip(edge.source+' → '+edge.target+'\nRepresentative changes: '+edge.count+'\nExact range: '+edge.minimum+'–'+edge.maximum)
        })));
        network.selectNodes(selected?[selected]:[],false);
        note.textContent='Node size: '+names[key]+'.'+(key==='source_hub_ratio'?' Sink ≈ 0 · hub = 0.5 · source ≈ 1.':'');
        if(selected){
          const connections=activeEdges.filter(edge=>edge.source===selected || edge.target===selected).length;
          const value=data.nodes.find(node=>node.id===selected).metrics[key];
          status.textContent=selected+': '+connections+' directed connection'+(connections===1?'':'s')+'; '+names[key].toLowerCase()+' '+format(value)+'. Direct neighbours are highlighted.';
        }else status.textContent=data.nodes.length+' countries · '+activeEdges.length+' directed links. Drag nodes, pan the background, or zoom to explore.';
        if(!data.nodes.length) status.textContent='No known countries are available in this view.';
      }
      listen(metric,'change',update);listen(country,'change',update);
      listen(possible,'change',()=>{update();network.stabilize(100);});
      for(const button of host.querySelectorAll('[data-network-action]')) listen(button,'click',()=>{
        const action=button.dataset.networkAction;
        if(action==='fit') fitNetwork();
        else if(action==='clear'){country.value='';update();}
        else network.moveTo({scale:Math.max(0.05,Math.min(8,network.getScale()*(action==='zoom-in'?1.25:0.8))),animation:false});
      });
      network.on('click',event=>{country.value=event.nodes.length?event.nodes[0]:'';update();});
      network.once('stabilized',()=>{network.stopSimulation();fitNetwork();});
      update();
    }
  }
  window.ChronoCladeNetworks={mount,destroy};
  if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',()=>mount());
  else mount();
})();
"""


def widget_assets():
    """Return CSS, offline library and initializer; append once to the report HTML."""
    library = files("chronoclade").joinpath("vendor", "vis-network.min.js").read_text()
    # Protect HTML script boundaries even if a future upstream bundle contains one.
    library = re.sub(r"</script", r"<\\/script", library, flags=re.IGNORECASE)
    return (
        "<style>"
        + _CSS
        + "</style>\n<script data-cc-network-library>"
        + library
        + "</script>\n<script>"
        + _INITIALIZER
        + "</script>"
    )
