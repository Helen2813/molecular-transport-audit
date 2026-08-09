import {
  CheckCircle2,
  FileCheck2,
  FolderLock,
  Users,
} from "lucide-react";

import type {
  AuditSummary,
  RunDetail,
} from "../../types/audit";

import StatCard from "./StatCard";


interface OverviewStatsProps {
  run: RunDetail;
  summary: AuditSummary | null;
}


function prettyStage(
  stage: string,
): string {
  return stage
    .replaceAll(
      "_",
      " ",
    )
    .replace(
      /\b\w/g,
      (letter) =>
        letter.toUpperCase(),
    );
}


export default function OverviewStats({
  run,
  summary,
}: OverviewStatsProps) {
  const passed =
    run.validated &&
    run.validation_status ===
      "PASS";

  const reference =
    summary?.run
      .reference_cohort ??
    run.reference_cohort ??
    "—";

  const external =
    summary?.run
      .external_cohort ??
    run.external_cohort ??
    "—";

  const outcomeLoaded =
    summary?.run
      .outcome_loaded ??
    run.outcome_loaded;

  const schemaVersion =
    summary?.schema_version ??
    run.schema_version;

  const statusValue =
    passed
      ? "PASS"
      : run.status.toUpperCase();

  let statusDescription =
    "Waiting to start";

  if (
    run.status === "running"
  ) {
    statusDescription =
      `Running: ${
        prettyStage(
          run.stage,
        )
      }`;
  }

  if (
    run.status === "failed"
  ) {
    statusDescription =
      "Scientific run failed";
  }

  if (passed) {
    statusDescription =
      "All required checks passed";
  }

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
        value={statusValue}
        description={
          statusDescription
        }
        icon={CheckCircle2}
        iconClassName={
          passed
            ? "text-emerald-400"
            : run.status ===
                "failed"
              ? "text-red-400"
              : "text-cyan-400"
        }
        iconBackgroundClassName={
          passed
            ? (
              "bg-emerald-400/[0.09] " +
              "ring-1 ring-emerald-400/10"
            )
            : run.status ===
                "failed"
              ? (
                "bg-red-400/[0.09] " +
                "ring-1 ring-red-400/10"
              )
              : (
                "bg-cyan-400/[0.09] " +
                "ring-1 ring-cyan-400/10"
              )
        }
        valueClassName={
          passed
            ? "text-emerald-400"
            : run.status ===
                "failed"
              ? "text-red-300"
              : "text-cyan-300"
        }
      />

      <StatCard
        label="Reference Cohort"
        value={reference}
        description="Canine reference atlas"
        icon={Users}
        iconClassName="text-blue-400"
        iconBackgroundClassName="bg-blue-500/[0.10] ring-1 ring-blue-400/10"
      />

      <StatCard
        label="External Cohort"
        value={external}
        description="External audit dataset"
        icon={Users}
        iconClassName="text-violet-400"
        iconBackgroundClassName="bg-violet-500/[0.10] ring-1 ring-violet-400/10"
      />

      <StatCard
        label="Outcome Loaded"
        value={
          outcomeLoaded === null
            ? "—"
            : outcomeLoaded
              ? "True"
              : "False"
        }
        description={
          outcomeLoaded
            ? "Outcome data available"
            : "No outcome data loaded"
        }
        icon={FolderLock}
        iconClassName="text-amber-400"
        iconBackgroundClassName="bg-amber-400/[0.09] ring-1 ring-amber-400/10"
        valueClassName={
          outcomeLoaded
            ? "text-white"
            : "text-amber-400"
        }
      />

      <StatCard
        label="Summary Contract"
        value={
          schemaVersion
            ? `v${schemaVersion}`
            : "Pending"
        }
        description={
          schemaVersion
            ? "Audit schema version"
            : (
              "Available after " +
              "verification"
            )
        }
        icon={FileCheck2}
        iconClassName="text-slate-300"
        iconBackgroundClassName="bg-slate-400/[0.08] ring-1 ring-slate-400/10"
      />
    </section>
  );
}
