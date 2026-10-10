"""Established offline report stylesheet and embedded font."""

import base64
from pathlib import Path


def embedded_font_css() -> str:
    """Embed the report display face so a saved report remains fully offline."""

    font = Path(__file__).resolve().parents[1] / "data" / "fonts" / "archivo-latin-variable.woff2"
    if not font.is_file():
        return ""
    encoded = base64.b64encode(font.read_bytes()).decode("ascii")
    return (
        "@font-face{font-family:'Archivo';font-style:normal;font-weight:100 900;"
        "font-display:swap;src:url(data:font/woff2;base64," + encoded + ") format('woff2');}"
    )


def report_styles() -> str:
    """Shared offline stylesheet for the existing ChronoClade report design."""
    return f"""
    {embedded_font_css()}
    :root {{ color-scheme:light; --ink:#17191f; --muted:#5f6470; --line:#c9ccd4; --paper:#fff; --wash:#eef0f4; --blue:#2855a6; --red:#d63c2f; --amber:#9a6500; --green:#14734f; --rail:248px; }}
    * {{ box-sizing:border-box; }}
    html {{ scroll-behavior:smooth; }}
    body {{ margin:0; background:var(--wash); color:var(--ink); font:16px/1.58 Archivo, "Helvetica Neue", sans-serif; }}
    ::selection {{ background:#ffdf78; color:var(--ink); }}
    a {{ color:var(--blue); text-underline-offset:3px; }} a:focus-visible,summary:focus-visible {{ outline:3px solid #ffb800; outline-offset:3px; }}
    .shell {{ width:min(1480px,100%); margin:auto; background:var(--paper); min-height:100vh; }}
    .identity {{ display:grid; grid-template-columns:var(--rail) minmax(0,1fr); border-bottom:2px solid var(--ink); }}
    .identity-mark {{ background:var(--ink); color:#fff; padding:30px 24px; font-weight:800; letter-spacing:.06em; }}
    .identity-copy {{ padding:24px 42px 26px; display:flex; flex-direction:column; gap:14px; align-items:flex-start; min-width:0; }}
    h1 {{ margin:0; font-size:clamp(32px,5vw,68px); line-height:1.12; letter-spacing:-.035em; max-width:900px; }} h1 i {{ font-style:italic; }} h1 span {{ display:block; }}
    .identity-copy p {{ margin:10px 0 0; color:var(--muted); max-width:70ch; }}
    .overall {{ width:max-content; max-width:100%; border:1px solid currentColor; padding:9px 12px; }} .overall.supported {{ color:var(--green); }} .overall.not_supported {{ color:var(--red); }} .overall small {{ color:currentColor; }}
    .overall small,.measure-strip small {{ display:block; text-transform:uppercase; letter-spacing:.08em; color:var(--muted); font-size:11px; font-weight:700; }}
    .overall b {{ display:block; font-size:21px; margin-top:4px; }}
    .contents {{ position:sticky; top:0; z-index:5; display:flex; gap:0; padding-left:var(--rail); background:#fff; border-bottom:1px solid var(--ink); overflow:auto; }}
    .contents a {{ flex:1; min-width:132px; padding:13px 16px; border-left:1px solid var(--line); color:var(--ink); text-decoration:none; font-size:13px; font-weight:700; }} .contents a:hover {{ background:#edf3ff; }}
    .stage {{ display:grid; grid-template-columns:var(--rail) minmax(0,1fr); border-bottom:2px solid var(--ink); scroll-margin-top:52px; }}
    .stage-index {{ padding:38px 24px; border-right:1px solid var(--ink); background:#f4f5f7; }}
    .stage-index span {{ display:block; font-size:74px; line-height:.8; font-weight:800; letter-spacing:-.06em; color:var(--blue); }}
    .stage-index b {{ display:block; margin-top:20px; text-transform:uppercase; letter-spacing:.1em; font-size:12px; }}
    .stage-body {{ padding:40px 44px 50px; min-width:0; }}
    .stage-head {{ display:flex; justify-content:space-between; gap:28px; align-items:start; margin-bottom:22px; }} .stage-head>div {{ min-width:0; }}
    h2 {{ margin:0; max-width:24ch; font-size:clamp(27px,3.2vw,44px); line-height:1.05; letter-spacing:-.025em; overflow-wrap:break-word; text-wrap:balance; }}
    h3 {{ margin:42px 0 8px; font-size:21px; }} p {{ max-width:75ch; }} .question {{ color:var(--muted); font-size:18px; margin:8px 0 0; }}
    .decision {{ display:block; min-width:150px; padding:10px 12px; border:2px solid currentColor; text-align:center; font-size:12px; letter-spacing:.08em; }}
    .decision.proceed {{ color:var(--green); }} .decision.stop {{ color:var(--red); }} .decision.review {{ color:var(--blue); }}
    figure {{ margin:28px 0 0; min-width:0; }} figure img {{ display:block; width:100%; height:auto; border:1px solid var(--ink); background:#fff; }} figcaption {{ margin-top:8px; color:var(--muted); font-size:13px; }}
    .evidence-layout {{ display:grid; grid-template-columns:minmax(0,1fr) 270px; gap:22px; align-items:start; min-width:0; }} .evidence-layout>* {{ min-width:0; }} .evidence-layout .evidence-files {{ grid-column:2; grid-row:1/4; margin-top:28px; min-width:0; }}
    .evidence-layout .download a {{ grid-template-columns:minmax(0,1fr) 42px; }} .evidence-layout .download small {{ grid-column:1/-1; grid-row:2; }}
    .figure-scroll {{ max-width:100%; overflow:auto; scrollbar-color:var(--blue) #e5e7ec; }} .figure-scroll:focus-visible,.table-scroll:focus-visible {{ outline:3px solid #ffb800; outline-offset:3px; }}
    .measure-strip {{ display:grid; grid-template-columns:repeat(3,1fr); border:1px solid var(--ink); margin:24px 0; }}
    .measure-strip > div {{ padding:16px; border-right:1px solid var(--line); }} .measure-strip > div:last-child {{ border:0; }}
    .measure-strip b {{ display:block; font-size:21px; margin:4px 0 2px; font-variant-numeric:tabular-nums; }} .measure-strip span {{ color:var(--muted); font-size:12px; }}
    .confidence-grid {{ display:grid; grid-template-columns:repeat(3,1fr); border:1px solid var(--ink); margin:24px 0; }} .confidence-grid>div {{ padding:16px; border-right:1px solid var(--line); border-bottom:1px solid var(--line); }} .confidence-grid>div:nth-child(3n) {{ border-right:0; }} .confidence-grid>div:nth-last-child(-n+3) {{ border-bottom:0; }} .confidence-grid small {{ display:block; text-transform:uppercase; letter-spacing:.08em; color:var(--muted); font-size:11px; font-weight:700; }} .confidence-grid b {{ display:block; font-size:21px; margin:4px 0 2px; font-variant-numeric:tabular-nums; }} .confidence-grid span {{ color:var(--muted); font-size:12px; }}
    .recombination-measures {{ grid-template-columns:repeat(5,minmax(0,1fr)); }} .recombination-measures>div {{ border-bottom:0; }} .recombination-measures>div:nth-child(3n) {{ border-right:1px solid var(--line); }} .recombination-measures>div:last-child {{ border-right:0; }}
    .logic {{ display:grid; grid-template-columns:1fr auto 1fr; align-items:center; gap:18px; margin:26px 0; }} .logic div {{ padding:20px; background:#f4f5f7; border:1px solid var(--line); }} .logic b {{ display:block; font-size:18px; }} .logic span {{ font-size:30px; color:var(--muted); }}
    details.evidence-files {{ margin-top:22px; border-top:1px solid var(--ink); border-bottom:1px solid var(--ink); }} details.evidence-files summary {{ padding:14px 0; cursor:pointer; font-weight:800; }}
    .downloads {{ list-style:none; margin:0 0 14px; padding:8px 12px 12px; border-top:1px solid var(--line); }} .download a {{ display:grid; grid-template-columns:minmax(180px,1fr) 2fr 54px; gap:14px; padding:11px 4px; border-bottom:1px solid var(--line); text-decoration:none; align-items:center; }} .download small {{ color:var(--muted); }} .download b {{ text-align:right; font-size:11px; letter-spacing:.08em; color:var(--muted); }}
    .scenario,.card,.verdict {{ margin:22px 0 0; padding:22px 0; border-top:1px solid var(--ink); background:#fff; }} .scenario strong {{ display:block; font-size:28px; max-width:34ch; }} .verdict.supported>strong {{ color:var(--green); }} .verdict.not_supported>strong {{ color:var(--red); }} .confidence,.evidence {{ font-weight:800; }} .guardrail,.caveat {{ background:#fff8dc; padding:15px; }} .card img {{ display:block; max-width:100%; height:auto; }} .check-list {{ list-style:none; margin:22px 0; padding:0; border-top:1px solid var(--line); }} .check-list li {{ display:grid; grid-template-columns:22px minmax(0,1fr); gap:10px; padding:12px 0; border-bottom:1px solid var(--line); }} .check-list li::before {{ content:"✓"; color:var(--green); font-weight:800; }}
    table {{ width:100%; border-collapse:collapse; margin-top:16px; font-size:14px; }} th,td {{ padding:10px 8px; text-align:left; border-bottom:1px solid var(--line); vertical-align:top; }} th {{ font-size:11px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }}
    .table-scroll {{ max-width:100%; overflow:auto; scrollbar-color:var(--blue) #e5e7ec; }} .scroll-hint {{ display:none; }}
    .report-close {{ display:grid; grid-template-columns:var(--rail) 1fr; background:var(--ink); color:#fff; }} .report-close b {{ padding:30px 24px; border-right:1px solid #555b67; }} .report-close div {{ padding:30px 42px; }} .report-close a {{ color:#fff; }}
    #root-to-tip .stage-body {{ display:grid; grid-template-columns:minmax(270px,.72fr) minmax(560px,1.55fr); column-gap:30px; align-items:start; }} #root-to-tip .stage-head,#root-to-tip .stage-body>p,#root-to-tip .measure-strip {{ grid-column:1; }} #root-to-tip .stage-head {{ display:block; }} #root-to-tip .decision {{ margin-top:18px; width:max-content; }} #root-to-tip .measure-strip {{ grid-template-columns:1fr; }} #root-to-tip .measure-strip>div {{ border-right:0; border-bottom:1px solid var(--line); }} #root-to-tip .evidence-layout {{ grid-column:2; grid-row:1/5; grid-template-columns:minmax(0,1fr) 210px; }} #root-to-tip .logic {{ grid-column:1; }}
    @media(max-width:1250px) {{ #root-to-tip .stage-body {{ display:block; }} .evidence-layout,#root-to-tip .evidence-layout {{ grid-template-columns:1fr; }} .evidence-layout .evidence-files {{ grid-column:1; grid-row:auto; margin-top:8px; }} }}
    @media(max-width:800px) {{ :root {{ --rail:76px; }} .identity-copy {{ display:block; padding:22px 20px; }} .overall {{ margin-top:22px; }} .contents {{ padding-left:0; overflow:visible; }} .contents a {{ min-width:0; padding:12px 7px; text-align:center; font-size:11px; white-space:nowrap; }} .stage-index {{ padding:28px 10px; }} .stage-index span {{ font-size:48px; }} .stage-index b {{ writing-mode:vertical-rl; margin:18px auto 0; }} .stage-body {{ padding:30px 18px 38px; }} .stage-head {{ display:block; }} .decision {{ margin-top:18px; width:max-content; max-width:100%; }} .measure-strip,.confidence-grid {{ grid-template-columns:1fr; }} .measure-strip > div,.confidence-grid>div {{ border-right:0; border-bottom:1px solid var(--line); min-width:0; }} .confidence-grid b {{ overflow-wrap:anywhere; }} .confidence-grid>div:nth-child(3n),.recombination-measures>div:nth-child(3n) {{ border-right:0; }} .confidence-grid>div:nth-last-child(-n+3) {{ border-bottom:1px solid var(--line); }} .confidence-grid>div:last-child {{ border-bottom:0; }} .logic {{ grid-template-columns:1fr; }} .logic > span {{ transform:rotate(90deg); justify-self:center; }} .download a {{ grid-template-columns:minmax(0,1fr) 42px; }} .download span,.download small {{ overflow-wrap:anywhere; }} .download small {{ grid-column:1/-1; grid-row:2; }} .scroll-hint {{ display:block; position:sticky; left:0; width:max-content; max-width:100%; margin-top:12px; padding:9px 10px; background:#edf3ff; color:var(--blue); font-size:12px; font-weight:700; }} .table-scroll:not(.compact-table) table {{ min-width:680px; }} .table-scroll:not(.compact-table) th:first-child,.table-scroll:not(.compact-table) td:first-child {{ position:sticky; left:0; background:#fff; z-index:1; }} .compact-table .scroll-hint,.compact-table thead {{ display:none; }} .compact-table table,.compact-table tbody,.compact-table tr,.compact-table td {{ display:block; min-width:0; width:100%; }} .compact-table tr {{ padding:10px 0; border-bottom:1px solid var(--line); }} .compact-table td {{ padding:3px 0; border:0; }} .compact-table td:first-child {{ font-weight:700; }} .figure-scroll img {{ min-width:760px; }} }}
    @media(prefers-reduced-motion:reduce) {{ html {{ scroll-behavior:auto; }} }}
    @media print {{ .contents {{ display:none; }} body,.shell {{ background:#fff; }} .stage {{ break-inside:avoid-page; }} details {{ display:block; }} details > * {{ display:block; }} .report-close {{ color:#000; background:#fff; border-top:2px solid #000; }} .report-close a {{ color:#000; }} }}
    .shell {{ width:min(1180px,100%); }}
    .identity {{ display:block; }} .identity-mark {{ padding:16px clamp(20px,5vw,54px); overflow-wrap:anywhere; }}
    .identity-copy {{ padding:30px clamp(20px,5vw,54px); }} h1 {{ font-size:clamp(32px,5vw,56px); }}
    .contents {{ padding-left:0; }} .contents a {{ flex:none; min-width:0; white-space:nowrap; }}
    .stage {{ display:block; }} .stage-index {{ display:none; }} .stage-body {{ padding:38px clamp(20px,5vw,54px); }}
    .stage .stage {{ border:0; }} .stage .stage .stage-body {{ padding:24px 0; }}
    #root-to-tip .stage-body {{ display:block; }} .evidence-layout,#root-to-tip .evidence-layout {{ display:block; }}
    .evidence-layout .evidence-files {{ margin-top:22px; }} .report-close {{ display:block; }}
    h2 {{ max-width:34ch; line-height:1.16; }} summary {{ cursor:pointer; padding:14px 0; font-weight:700; }}
    .context-geography {{ min-width:0; max-width:100%; }} .context-geography img,.context-geography svg {{ width:100%; height:auto; }}
    .demonstration strong {{ display:block; }} .demonstration p {{ margin:8px 0; }}
    @media(max-width:800px) {{
      .neighbourhood [aria-label="Closest relatives by final clonal SNPs"] table {{ min-width:0; }}
      .neighbourhood [aria-label="Closest relatives by final clonal SNPs"] thead {{ display:none; }}
      .neighbourhood [aria-label="Closest relatives by final clonal SNPs"] tr {{ display:block; padding:12px 0; border-bottom:1px solid var(--line); }}
      .neighbourhood [aria-label="Closest relatives by final clonal SNPs"] td {{ display:grid; grid-template-columns:140px minmax(0,1fr); gap:12px; padding:5px 0; border:0; overflow-wrap:anywhere; position:static; }}
      .neighbourhood [aria-label="Closest relatives by final clonal SNPs"] td::before {{ content:attr(data-label); font-weight:700; }}
    }}
    @media(max-width:800px) {{ .contents {{ overflow:auto; }} .contents a {{ flex:none; padding:12px 14px; font-size:13px; }} .stage-body {{ padding:30px 20px; }} .identity-copy {{ padding:28px 20px; }} }}
    """
