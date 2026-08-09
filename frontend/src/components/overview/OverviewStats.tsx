import {
  CheckCircle2,
  FileCheck2,
  FolderLock,
  Users,
} from "lucide-react";

import StatCard from "./StatCard";

import type {
  AuditSummary,
  RunDetail,
} from "../../types/audit";

interface OverviewStatsProps {
  run: RunDetail;
  summary: AuditSummary;
}

export default function OverviewStats({
  run,
  summary,
}: OverviewStatsProps) {
  const passed =
    run.validated &&
    run.validation_status === "PASS";

  return (
    <section
      className="
        grid
        grid-cols-1
        gap-3
        sm:grid-cols-2
        xl:grid-cols-[0.95fr_1.1fr_1.05fr_0.95fr_0.9fr]
      "
    >
      <StatCard
        label="Status"
        value={
          passed
            ? "PASS"
            : run.status.toUpperCase()
        }
        description={
          passed
            ? "All required checks passed"
            : `Current stage: ${run.stage}`
        }
        icon={CheckCircle2}
        iconClassName={
          passed
            ? "text-emerald-400"
            : "text-amber-400"
        }
        iconBackgroundClassName={
          passed
            ? "bg-emerald-400/[0.09] ring-1 ring-emerald-400/10"
            : "bg-amber-400/[0.09] ring-1 ring-amber-400/10"
        }
        valueClassName={
          passed
            ? "text-emerald-400"
            : "text-amber-300"
        }
      />

      <StatCard
        label="Reference Cohort"
        value={
          summary.run.reference_cohort
        }
        description="Canine reference atlas"
        icon={Users}
        iconClassName="text-blue-400"
        iconBackgroundClassName="bg-blue-500/[0.10] ring-1 ring-blue-400/10"
      />

      <StatCard
        label="External Cohort"
        value={
          summary.run.external_cohort
        }
        description="External audit dataset"
        icon={Users}
        iconClassName="text-violet-400"
        iconBackgroundClassName="bg-violet-500/[0.10] ring-1 ring-violet-400/10"
      />

      <StatCard
        label="Outcome Loaded"
        value={
          summary.run.outcome_loaded
            ? "True"
            : "False"
        }
        description={
          summary.run.outcome_loaded
            ? "Outcome data available"
            : "No outcome data loaded"
        }
        icon={FolderLock}
        iconClassName="text-amber-400"
        iconBackgroundClassName="bg-amber-400/[0.09] ring-1 ring-amber-400/10"
        valueClassName={
          summary.run.outcome_loaded
            ? "text-white"
            : "text-amber-400"
        }
      />

      <StatCard
        label="Summary Contract"
        value={`v${summary.schema_version}`}
        description="Audit schema version"
        icon={FileCheck2}
        iconClassName="text-slate-300"
        iconBackgroundClassName="bg-slate-400/[0.08] ring-1 ring-slate-400/10"
      />
    </section>
  );
}
