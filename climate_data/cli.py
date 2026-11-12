"""资料摘要、分位增量订正与锁定模板秩恢复的命令入口。

所有命令只读取来源文件，绝不回写；inspect 输出口径与历史版本一致。
"""
import argparse
import json
import sys

from .io import load
from .domain import summarize
from .correction import correct_document, apply_restoration, TemplateRestorationError


def _emit(payload):
    print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="合成降水情景的分位增量订正与同期秩恢复")
    parser.add_argument('command', choices=['inspect', 'correct', 'restore'])
    parser.add_argument('file')
    args = parser.parse_args(argv)

    value, digest = load(args.file)

    if args.command == 'inspect':
        _emit({'source_sha256': digest, 'summary': summarize(value)})
        return 0

    corrected = correct_document(value)
    if args.command == 'correct':
        _emit({'source_sha256': digest, 'correction': corrected})
        return 0

    try:
        restored = apply_restoration(value, corrected)
    except TemplateRestorationError as exc:
        # 模板不合规时整次恢复拒绝：仍给出订正版，明确不执行恢复。
        _emit({
            'source_sha256': digest,
            'restoration': {
                'accepted': False,
                'reason': str(exc),
                'correction': corrected,
            },
        })
        return 3
    _emit({'source_sha256': digest,
           'correction': corrected,
           'restoration': restored})
    return 0


if __name__ == '__main__':
    sys.exit(main())
