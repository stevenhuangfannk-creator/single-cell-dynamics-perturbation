"""Build a GitHub gallery from actual module figure manifests."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def build():
    records = []
    for module in sorted(ROOT.glob('modules/*')):
        path = module / 'figures/figure_manifest.json'
        if not path.exists():
            path = module / 'results/figure_manifest.json'
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding='utf-8-sig'))
        records.extend(payload['figures'] if isinstance(payload, dict) else payload)
    lines = ['# 真实图件索引 / Figure gallery', '',
             '每张预览均来自保存的计算结果。图中距离、时间或网络倍数须按对应方法说明解释；没有实际结果的面板不以示意图填充。单图重绘只读取结果对象，不重新训练模型。', '',
             '| 编号、方法与代码入口 | 数据、科学含义与状态 | 真实图件预览及导出 |',
             '|---|---|---|']
    for item in records:
        required = ('id', 'title_cn', 'meaning', 'dataset', 'reference', 'analysis_script', 'plot_script', 'png', 'pdf', 'svg', 'status', 'validity')
        missing = [key for key in required if key not in item]
        if missing:
            raise ValueError(f'{item.get("id")}: manifest missing {missing}')
        for key in ('analysis_script', 'plot_script', 'png', 'pdf', 'svg'):
            path = ROOT / (item[key].split(' ', 1)[0] if key.endswith('script') else item[key])
            if not path.is_file() or not path.stat().st_size:
                raise FileNotFoundError(f'{item["id"]}: {key} {path}')
        analysis_scripts = item.get('analysis_scripts', [item['analysis_script']])
        for script in analysis_scripts:
            if not (ROOT / script).is_file():
                raise FileNotFoundError(f'{item["id"]}: analysis script {script}')
        def clean(value):
            return str(value).replace('|', '\\|').replace('\n', '<br>')
        analysis = item['analysis_script'].split(' ', 1)[0]
        plotting = item['plot_script'].split(' ', 1)[0]
        method = item.get('method', 'scTenifoldKnk' if '01_scTenifoldKnk' in analysis else 'scTour')
        command = item.get('redraw_command', item['plot_script'])
        analysis_links = ' · '.join(f'[分析{index if len(analysis_scripts) > 1 else ""}](../{script})'
                                    for index, script in enumerate(analysis_scripts, 1))
        left = f'**{clean(item["id"])} · {clean(item["title_cn"])}**<br>{clean(method)}<br>{analysis_links} · [重绘](../{plotting})<br><code>{clean(command)}</code><br>[参考来源]({item["reference"]})'
        middle = f'{clean(item["dataset"])}<br>{clean(item["meaning"])}<br>**{item["status"]}**：{clean(item["validity"])}'
        right = f'<img src="../{item["png"]}" width="260" alt="{clean(item["id"])}"><br>[PNG](../{item["png"]}) · [PDF](../{item["pdf"]}) · [SVG](../{item["svg"]})'
        lines.append(f'| {left} | {middle} | {right} |')
    output = ROOT / 'gallery/FIGURE_INDEX_CN.md'
    output.parent.mkdir(exist_ok=True)
    output.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Indexed {len(records)} actual figures: {output.relative_to(ROOT)}')
    return records


if __name__ == '__main__':
    build()
