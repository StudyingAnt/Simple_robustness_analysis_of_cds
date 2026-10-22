"""Summarize a completed pilot without generating any additional trajectories."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, help='Earlier pilot for prefix and precision comparison')
    args = parser.parse_args()
    directory = args.results
    manifest = json.loads((directory / 'manifest.json').read_text())
    rows = json.loads((directory / 'summary.json').read_text())
    n = manifest['n_traj_per_case']
    assert isinstance(n, int) and n > 0
    assert manifest['repeats'] == 1 and manifest['trajectory_warmup_count'] == 0
    assert manifest['successful_cases'] == 36 and manifest['failed_cases'] == 0
    assert len(rows) == 36
    expected = {(item['transcript_id'], signature)
                for item in manifest['selection'] for signature in manifest['signatures']}
    assert {(row['transcript_id'], row['signature']) for row in rows} == expected
    audited = 0
    total_steps = 0
    for row in rows:
        case = json.loads((directory / 'cases' / f"{row['transcript_id']}_{row['signature']}.json").read_text())
        assert len(case['paths']) == n and case['all_path_checks_passed']
        assert case['n_complete'] == n and case['n_incomplete'] == 0
        assert [path['trajectory_id'] for path in case['paths']] == list(range(n))
        assert all(path['terminal_reason'] == 'nonsynonymous' for path in case['paths'])
        accepted = [path['accepted_steps'] for path in case['paths']]
        assert statistics.mean(accepted) == row['mean_accepted']
        assert all(path['n_steps'] == path['accepted_steps'] + 1 for path in case['paths'])
        audited += len(case['paths'])
        total_steps += sum(path['n_steps'] for path in case['paths'])
    assert audited == 36 * n
    comparison = None
    if args.baseline:
        before_manifest = json.loads((args.baseline / 'manifest.json').read_text())
        before_rows = json.loads((args.baseline / 'summary.json').read_text())
        before_n = before_manifest['n_traj_per_case']
        assert before_n < n
        assert before_manifest['selected_fasta_sha256'] == manifest['selected_fasta_sha256']
        assert before_manifest['package_commit'] == manifest['package_commit']
        assert before_manifest['signatures'] == manifest['signatures']
        prefix_count = 0
        differences, se_ratios = [], []
        zero_baseline_now_positive = []
        for row in rows:
            name = f"{row['transcript_id']}_{row['signature']}.json"
            before = json.loads((args.baseline / 'cases' / name).read_text())
            after = json.loads((directory / 'cases' / name).read_text())
            assert before['case_seed'] == after['case_seed']
            assert before['paths'] == after['paths'][:before_n]
            prefix_count += before_n
            differences.append(abs(after['mean_accepted'] - before['mean_accepted']))
            if before['accepted_mc_se'] > 0:
                se_ratios.append(after['accepted_mc_se'] / before['accepted_mc_se'])
            elif after['mean_accepted'] > 0:
                zero_baseline_now_positive.append(dict(gene=row['gene'], signature=row['signature'],
                                                       before=before['mean_accepted'], after=after['mean_accepted']))
        comparison = dict(
            baseline_n=before_n, current_n=n, matching_prefix_paths=prefix_count,
            all_prefix_paths_identical=True,
            baseline_total_driver_wall_seconds=before_manifest['total_driver_wall_seconds'],
            baseline_simulation_seconds=sum(row['simulation_seconds'] for row in before_rows),
            baseline_median_mc_se=statistics.median(row['accepted_mc_se'] for row in before_rows),
            current_median_mc_se=statistics.median(row['accepted_mc_se'] for row in rows),
            median_mc_se_ratio_for_positive_baseline=statistics.median(se_ratios),
            median_absolute_mean_change=statistics.median(differences),
            maximum_absolute_mean_change=max(differences),
            baseline_zero_mean_now_positive=zero_baseline_now_positive,
            means_are_nested_samples=True,
        )
    grouped = []
    for name in ['short', 'medium', 'long', 'very_long']:
        group = [row for row in rows if row['bin'] == name]
        grouped.append(dict(bin=name, population=manifest['population_by_bin'][name],
                            sampled_transcripts=len({row['transcript_id'] for row in group}),
                            sampled_cds_min=min(row['cds_nt'] for row in group),
                            sampled_cds_max=max(row['cds_nt'] for row in group),
                            seconds_min=min(row['simulation_seconds'] for row in group),
                            seconds_median=statistics.median(row['simulation_seconds'] for row in group),
                            seconds_max=max(row['simulation_seconds'] for row in group),
                            peak_worker_rss_mib_max=max(row['peak_rss_mib'] for row in group)))
    simulation = sum(row['simulation_seconds'] for row in rows)
    validation = sum(row['validation_seconds'] for row in rows)
    ttn_seconds = sum(row['simulation_seconds'] for row in rows if row['gene'] == 'TTN')
    projection_target = 1000 if n >= 100 else 100
    assessment = dict(
        measured_n_traj_per_case=n, measured_total_trajectories=audited,
        total_mutation_steps=total_steps, all_trajectory_audits_passed=True,
        simulation_seconds=simulation,
        simulation_cpu_seconds=sum(row['simulation_cpu_seconds'] for row in rows),
        validation_seconds=validation, total_driver_wall_seconds=manifest['total_driver_wall_seconds'],
        worker_rss_mib_min=min(row['peak_rss_mib'] for row in rows),
        worker_rss_mib_max=max(row['peak_rss_mib'] for row in rows),
        ttn_fraction_of_sample_simulation_time=ttn_seconds / simulation,
        median_accepted_mc_se=statistics.median(row['accepted_mc_se'] for row in rows),
        maximum_accepted_mc_se=max(row['accepted_mc_se'] for row in rows),
        projected_n_traj_per_case=projection_target,
        estimated_larger_pilot_seconds=manifest['total_driver_wall_seconds'] +
            (projection_target / n - 1) * (simulation + validation),
        projection_method='Fixed overhead held constant; simulation and validation scaled by trajectory count.',
        projection_is_measured=False, repeated_timing_runs_per_case=1,
        comparison_with_baseline=comparison,
        groups=grouped,
        cpu_model=next(line.split(':', 1)[1].strip()
                       for line in Path('/proc/cpuinfo').read_text().splitlines()
                       if line.startswith('model name')),
    )
    (directory / 'assessment.json').write_text(json.dumps(assessment, indent=2) + '\n')
    with (directory / 'length_groups.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(grouped[0]))
        writer.writeheader()
        writer.writerows(grouped)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
    colors = dict(SBS1='#3366aa', SBS2='#d05a28', SBS13='#218563')
    for signature, color in colors.items():
        group = sorted((row for row in rows if row['signature'] == signature),
                       key=lambda row: row['cds_nt'])
        lengths = [row['cds_nt'] for row in group]
        axes[0].plot(lengths, [row['simulation_seconds'] for row in group],
                     marker='o', markersize=4, linewidth=1, label=signature, color=color)
        axes[1].plot(lengths, [row['peak_rss_mib'] for row in group],
                     marker='o', markersize=4, linewidth=1, label=signature, color=color)
    for axis in axes:
        axis.set_xscale('log')
        axis.set_xlabel('CDS length (nt)')
        axis.grid(True, which='major', alpha=.25)
    axes[0].set_yscale('log')
    axes[0].set_ylabel(f'Simulation time for {n} trajectories (s)')
    axes[0].set_title('Time: profile cache prepared before timing')
    axes[0].legend(frameon=False)
    axes[1].set_ylabel('Peak worker RSS (MiB)')
    axes[1].set_title('Memory: imports, simulation and validation')
    fig.suptitle(f'Sinabro 3.0.0 pilot: 12 transcripts, 3 signatures, {n} trajectories each')
    fig.savefig(directory / 'length_time_memory.png', dpi=200)
    fig.savefig(directory / 'length_time_memory.svg')
    plt.close(fig)

    lines = [
        f'# {n}개 경로 파일럿 결과',
        '',
        f'12개 전사체 × SBS1/SBS2/SBS13 × 조합당 {n}개 = {audited:,}개 경로를 실행했다. '
        f'36개 조합 모두 정상 완료했고, {audited:,}개 모두 첫 단백질 변화로 종료했다. '
        '추가 난수 반복은 수행하지 않았다.',
        '',
        '논문 원본 simulate_single_gene.py의 n_traj는 1,000이다. '
        f'이번에는 조합당 반복 경로 수를 {n}개로 낮췄다. '
        '각 경로는 첫 아미노산 변화까지 진행하며, RUMS는 마지막 비동의적 변이를 제외한 변이 수의 평균이다.',
        '',
        f"실제 시뮬레이션 합계: {simulation:.3f}초. "
        f"표본 선정·각 작업 프로세스 시작·검증·저장을 포함한 전체: {manifest['total_driver_wall_seconds']:.3f}초.",
        '',
        f'| 구간 | 전체 전사체 수 | 표본 CDS 길이 (nt) | 조합당 {n}개 계산 시간 범위 (초) |',
        '|---|---:|---:|---:|',
    ]
    for group in grouped:
        lines.append(f"| {group['bin']} | {group['population']:,} | "
                     f"{group['sampled_cds_min']:,}–{group['sampled_cds_max']:,} | "
                     f"{group['seconds_min']:.3f}–{group['seconds_max']:.3f} |")
    lines += [
        '',
        f"작업 프로세스의 최고 RSS는 {assessment['worker_rss_mib_min']:.1f}–"
        f"{assessment['worker_rss_mib_max']:.1f} MiB였다. 이 수치는 각 worker의 import·계산·검증 중 최고값이다.",
        '',
        '![길이별 시간과 메모리](length_time_memory.png)',
        '',
        '## 평가',
        '',
        '- 이번 범위에서는 입력·프레임·단일 치환·종료 조건 오류나 미완료 경로를 발견하지 못했다. '
        '각 경로의 모든 중간 단백질, 마지막 변화, cDNA 좌표, flanking base와 accepted count를 검증했다.',
        '- 계산 시간은 대체로 CDS 길이와 실제 누적 변이 단계 수에 따라 증가한다. '
        '가장 긴 TTN은 107,976 nt다.',
        f"- TTN 세 조합이 표본 계산 시간의 {100 * ttn_seconds / simulation:.1f}%를 차지했다. "
        '길이 구간을 같은 개수로 뽑고 최장 서열도 포함했으므로 이 평균을 전체 전사체 수에 그대로 곱하면 안 된다.',
        f"- RUMS 평균의 Monte Carlo 표준오차는 조합별 중앙값 {assessment['median_accepted_mc_se']:.3f}, "
        f"최대 {assessment['maximum_accepted_mc_se']:.3f}였다. "
        '이는 Monte Carlo 표준오차이며 시간 측정의 오차가 아니다.',
        f"- 동일 36개 조합을 {projection_target:,}개씩 돌린다면 단순 선형 환산은 약 "
        f"{assessment['estimated_larger_pilot_seconds']:.0f}초다. "
        '고정 비용을 유지한 추정이며 실제 측정값이나 신뢰구간이 아니다. 전체 유전체 실행의 추정도 아니다.',
        '',
        '## 조건과 한계',
        '',
        f"- Python {manifest['python']}, NumPy {manifest['numpy']}, Biopython {manifest['biopython']}; Sinabro 3.0.0.",
        '- COSMIC v3.3.1 GRCh37, protein-identity 기준, CDS 시작 1, 고정 master seed 42.',
        '- 한 번에 worker 1개, 수치 라이브러리 스레드 1개. 각 조합은 새 프로세스에서 실행했다.',
        '- 표본은 1,000 / 3,000 / 10,000 nt 경계로 나누고 서로 다른 유전자를 선정했다. '
        '앞 세 구간은 길이의 10/50/90 분위수 근처, 마지막 구간은 10/50 분위수와 최장을 택했다.',
        '- 내부 표준코드 종결 코돈을 가진 73개는 이번 timing 표본에서만 제외했다. 원본 전처리 결과는 유지했다.',
        '- 프로파일 로딩을 계산 전에 준비했으며 그 비용은 setup과 전체 시간에 포함된다. '
        '정확성 검증 시간은 simulation 시간에서 분리했다.',
        '- 각 조합의 시간 반복 측정은 1회다. 기계 부하에 따른 시간 변동과 병렬 확장성은 아직 평가하지 않았다.',
        '- 세 시그니처와 strict 기준의 결과다. 다른 76개 SBS 시그니처나 BLOSUM 기준의 시간은 이 결과로 확정하지 않는다.',
    ]
    if comparison:
        lines += [
            '',
            f"이전 {comparison['baseline_n']}개 측정과 비교:",
            '',
            f"- 동일 조합의 앞 {comparison['baseline_n']}개 경로, 총 "
            f"{comparison['matching_prefix_paths']:,}개가 이전과 모두 동일하다. "
            '단계 수, 종료 이유, accepted count 및 전체 DNA 경로 해시를 대조했다.',
            f"- Monte Carlo 표준오차의 조합별 중앙값: {comparison['baseline_median_mc_se']:.3f} → "
            f"{comparison['current_median_mc_se']:.3f}.",
            f"- RUMS 평균 변화의 절댓값: 중앙값 {comparison['median_absolute_mean_change']:.3f}, "
            f"최대 {comparison['maximum_absolute_mean_change']:.3f}.",
            '- 두 표본은 앞부분을 공유한다. 독립적인 반복 실험으로 해석하지 않는다.',
            '',
        ]
        if comparison['baseline_zero_mean_now_positive']:
            lines.append('이전에는 0이었으나 추가 경로에서 양수가 관찰된 조합: ' +
                         json.dumps(comparison['baseline_zero_mean_now_positive'], ensure_ascii=False))
            lines.append('')
    lines += [
        '',
        '## 세부 결과',
        '',
        '| Gene | CDS nt | SBS1 초 / RUMS | SBS2 초 / RUMS | SBS13 초 / RUMS |',
        '|---|---:|---:|---:|---:|',
    ]
    for selected in manifest['selection']:
        cells = []
        for signature in manifest['signatures']:
            row = next(row for row in rows if row['transcript_id'] == selected['transcript_id']
                       and row['signature'] == signature)
            cells.append(f"{row['simulation_seconds']:.3f} / {row['mean_accepted']:.2f}")
        lines.append(f"| {selected['gene']} | {selected['cds_nt']:,} | " + ' | '.join(cells) + ' |')
    lines += ['', '원시 경로 수·종료 상태·seed·해시: cases/의 조합별 JSON. '
              '전체 수치 표: summary.csv. 환경·입력·소스 기록: manifest.json. '
              '이 평가 스크립트는 추가 시뮬레이션을 하지 않는다.', '']
    (directory / 'report.md').write_text('\n'.join(lines))
    print(json.dumps(assessment, indent=2))


if __name__ == '__main__':
    main()
