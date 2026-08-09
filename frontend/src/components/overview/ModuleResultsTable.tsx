import {
  Download,
  Info,
} from "lucide-react";

import type {
  ModuleResult,
} from "../../types/audit";


interface ModuleResultsTableProps {
  modules: ModuleResult[];
}


function formatMetric(
  value: number | null,
): string {
  if (value === null) {
    return "—";
  }

  return value.toFixed(6);
}


function formatEmpiricalP(
  value: number | null,
): string {
  if (value === null) {
    return "—";
  }

  return value.toFixed(6);
}


function classificationLabel(
  value: ModuleResult[
    "classification"
  ]["short"],
): string {
  const labels = {
    strong: "STRONG",
    partial: "PARTIAL",
    limited: "LIMITED",
    no_clear: "NO_CLEAR",
  };

  return labels[value];
}


function ClassificationBadge({
  value,
}: {
  value: ModuleResult[
    "classification"
  ]["short"];
}) {
  const styles = {
    strong:
      "border-emerald-400/15 bg-emerald-400/[0.07] text-emerald-400",
    partial:
      "border-blue-400/15 bg-blue-400/[0.07] text-blue-300",
    limited:
      "border-violet-400/15 bg-violet-400/[0.07] text-violet-300",
    no_clear:
      "border-amber-400/10 bg-amber-400/[0.06] text-slate-300",
  };

  const dotStyles = {
    strong:
      "bg-emerald-400 shadow-[0_0_8px_rgba(52,211,153,0.65)]",
    partial:
      "bg-blue-400",
    limited:
      "bg-violet-400",
    no_clear:
      "bg-amber-400",
  };

  return (
    <span
      className={[
        "inline-flex h-7 items-center gap-2 rounded-md border px-2.5 text-[11px] font-semibold tracking-[0.02em]",
        styles[value],
      ].join(" ")}
    >
      <span
        className={[
          "h-1.5 w-1.5 rounded-full",
          dotStyles[value],
        ].join(" ")}
      />

      {classificationLabel(value)}
    </span>
  );
}


function empiricalPClass(
  value: number | null,
): string {
  if (value === null) {
    return "text-slate-500";
  }

  if (value <= 0.01) {
    return "text-emerald-400";
  }

  return "text-amber-400";
}


function escapeCsv(
  value: string | number | null,
): string {
  if (value === null) {
    return "";
  }

  const text = String(value);

  if (
    text.includes(",") ||
    text.includes('"') ||
    text.includes("\n")
  ) {
    return `"${text.replaceAll(
      '"',
      '""',
    )}"`;
  }

  return text;
}


function exportModulesCsv(
  modules: ModuleResult[],
) {
  const header = [
    "module",
    "classification",
    "edge_spearman",
    "loading_spearman",
    "split_half_median",
    "random_panel_empirical_p",
  ];

  const rows = modules.map(
    (module) => [
      module.module,
      module.classification.short,
      module.direct_preservation
        .edge_spearman,
      module.direct_preservation
        .loading_spearman,
      module.reliability
        .split_half_median,
      module.random_control
        .empirical_p,
    ],
  );

  const csv = [
    header,
    ...rows,
  ]
    .map((row) =>
      row
        .map(escapeCsv)
        .join(","),
    )
    .join("\n");

  const blob = new Blob(
    [csv],
    {
      type:
        "text/csv;charset=utf-8",
    },
  );

  const url =
    URL.createObjectURL(blob);

  const link =
    document.createElement("a");

  link.href = url;
  link.download =
    "module_summary.csv";

  document.body.appendChild(
    link,
  );

  link.click();
  link.remove();

  URL.revokeObjectURL(url);
}


export default function ModuleResultsTable({
  modules,
}: ModuleResultsTableProps) {
  return (
    <section className="overflow-hidden rounded-[14px] border border-white/[0.07] bg-[#091725]/90 shadow-[0_12px_40px_rgba(0,0,0,0.12)]">
      <div className="flex h-[52px] items-center justify-between border-b border-white/[0.07] px-5">
        <div className="flex items-center gap-2">
          <h2 className="text-[14px] font-semibold tracking-[-0.015em] text-slate-100">
            Module Results
          </h2>

          <Info
            size={14}
            strokeWidth={1.8}
            className="text-slate-500"
          />
        </div>

        <button
          type="button"
          onClick={() =>
            exportModulesCsv(
              modules,
            )
          }
          className="flex h-8 items-center gap-2 rounded-lg border border-white/[0.06] bg-white/[0.025] px-3 text-[12px] font-medium text-slate-300 transition hover:border-white/[0.11] hover:bg-white/[0.05] hover:text-white"
        >
          <Download
            size={14}
            strokeWidth={1.8}
          />

          Export CSV
        </button>
      </div>

      <div className="overflow-x-auto">
        <table className="w-full min-w-[900px] border-collapse">
          <thead>
            <tr className="border-b border-white/[0.06] bg-white/[0.008]">
              <th className="w-[10%] px-5 py-3 text-left text-[12px] font-medium text-slate-400">
                Module
              </th>

              <th className="w-[17%] px-4 py-3 text-left text-[12px] font-medium text-slate-400">
                Classification
              </th>

              <th className="w-[18%] px-4 py-3 text-right text-[12px] font-medium text-slate-400">
                Edge Spearman ρ
              </th>

              <th className="w-[18%] px-4 py-3 text-right text-[12px] font-medium text-slate-400">
                Loading Spearman ρ
              </th>

              <th className="w-[18%] px-4 py-3 text-right text-[12px] font-medium text-slate-400">
                Split-Half Median ρ
              </th>

              <th className="w-[19%] px-5 py-3 text-right text-[12px] font-medium text-slate-400">
                Random Panel Empirical p
              </th>
            </tr>
          </thead>

          <tbody>
            {modules.map(
              (module) => (
                <tr
                  key={
                    module.module
                  }
                  className="group border-b border-white/[0.045] transition last:border-b-0 hover:bg-white/[0.018]"
                >
                  <td className="px-5 py-[12px]">
                    <span className="text-[14px] font-semibold text-slate-100">
                      {module.module}
                    </span>
                  </td>

                  <td className="px-4 py-[12px]">
                    <ClassificationBadge
                      value={
                        module
                          .classification
                          .short
                      }
                    />
                  </td>

                  <td className="px-4 py-[12px] text-right font-mono text-[14px] tabular-nums text-slate-200">
                    {formatMetric(
                      module
                        .direct_preservation
                        .edge_spearman,
                    )}
                  </td>

                  <td className="px-4 py-[12px] text-right font-mono text-[14px] tabular-nums text-slate-200">
                    {formatMetric(
                      module
                        .direct_preservation
                        .loading_spearman,
                    )}
                  </td>

                  <td className="px-4 py-[12px] text-right font-mono text-[14px] tabular-nums text-slate-200">
                    {formatMetric(
                      module
                        .reliability
                        .split_half_median,
                    )}
                  </td>

                  <td
                    className={[
                      "px-5 py-[12px] text-right font-mono text-[14px] font-medium tabular-nums",
                      empiricalPClass(
                        module
                          .random_control
                          .empirical_p,
                      ),
                    ].join(" ")}
                  >
                    {formatEmpiricalP(
                      module
                        .random_control
                        .empirical_p,
                    )}
                  </td>
                </tr>
              ),
            )}
          </tbody>
        </table>
      </div>
    </section>
  );
}
