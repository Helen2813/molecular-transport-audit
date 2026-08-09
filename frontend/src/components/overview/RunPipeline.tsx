import {
  BadgeCheck,
  Check,
  Circle,
  Clock3,
  FileCheck2,
  LoaderCircle,
  Network,
  ShieldCheck,
  Shuffle,
  Sigma,
  XCircle,
  type LucideIcon,
} from "lucide-react";

import {
  useEffect,
  useMemo,
  useState,
} from "react";

import type {
  RunDetail,
  RunStage,
} from "../../types/audit";


interface RunPipelineProps {
  run: RunDetail;
}


interface PipelineStep {
  stage: RunStage;
  label: string;
  description: string;
  icon: LucideIcon;
}


const PIPELINE_STEPS: PipelineStep[] = [
  {
    stage: "direct_preservation",
    label: "Direct Preservation",
    description: "Structure & coverage",
    icon: Network,
  },
  {
    stage: "permutation_inference",
    label: "Permutation Inference",
    description: "Statistical inference",
    icon: Sigma,
  },
  {
    stage: "reliability",
    label: "Reliability",
    description: "Split-half & LOO",
    icon: ShieldCheck,
  },
  {
    stage: "random_controls",
    label: "Random Controls",
    description: "Matched specificity",
    icon: Shuffle,
  },
  {
    stage: "summary_assembly",
    label: "Summary Assembly",
    description: "Contract generation",
    icon: FileCheck2,
  },
  {
    stage: "verification",
    label: "Verification",
    description: "Independent checks",
    icon: BadgeCheck,
  },
];


type VisualState =
  | "complete"
  | "active"
  | "pending";


function useCurrentTime(
  enabled: boolean,
): number {
  const [
    currentTime,
    setCurrentTime,
  ] = useState(
    Date.now(),
  );

  useEffect(() => {
    if (!enabled) {
      return;
    }

    const interval =
      window.setInterval(
        () => {
          setCurrentTime(
            Date.now(),
          );
        },
        1000,
      );

    return () => {
      window.clearInterval(
        interval,
      );
    };
  }, [enabled]);

  return currentTime;
}


function formatDuration(
  startedAt: string | null,
  finishedAt: string | null,
  now: number,
): string {
  if (!startedAt) {
    return "—";
  }

  const start =
    new Date(
      startedAt,
    ).getTime();

  const end = finishedAt
    ? new Date(
        finishedAt,
      ).getTime()
    : now;

  if (
    !Number.isFinite(start) ||
    !Number.isFinite(end)
  ) {
    return "—";
  }

  const seconds = Math.max(
    0,
    Math.floor(
      (end - start) / 1000,
    ),
  );

  const hours = Math.floor(
    seconds / 3600,
  );

  const minutes = Math.floor(
    (seconds % 3600) / 60,
  );

  const remainingSeconds =
    seconds % 60;

  return [
    hours,
    minutes,
    remainingSeconds,
  ]
    .map((value) =>
      String(value).padStart(
        2,
        "0",
      ),
    )
    .join(":");
}


function statusLabel(
  run: RunDetail,
): string {
  if (
    run.status === "completed" &&
    run.validated
  ) {
    return "Completed";
  }

  if (run.status === "running") {
    return "Running";
  }

  if (run.status === "queued") {
    return "Queued";
  }

  return "Failed";
}


function statusClasses(
  run: RunDetail,
): string {
  if (
    run.status === "completed" &&
    run.validated
  ) {
    return (
      "border-emerald-400/15 " +
      "bg-emerald-400/[0.07] " +
      "text-emerald-400"
    );
  }

  if (run.status === "running") {
    return (
      "border-cyan-400/15 " +
      "bg-cyan-400/[0.07] " +
      "text-cyan-300"
    );
  }

  if (run.status === "failed") {
    return (
      "border-red-400/15 " +
      "bg-red-400/[0.07] " +
      "text-red-300"
    );
  }

  return (
    "border-amber-400/15 " +
    "bg-amber-400/[0.07] " +
    "text-amber-300"
  );
}


function getVisualState(
  run: RunDetail,
  stepIndex: number,
  currentIndex: number,
): VisualState {
  if (
    run.status === "completed"
  ) {
    return "complete";
  }

  if (
    run.status !== "running"
  ) {
    return "pending";
  }

  if (currentIndex < 0) {
    return "pending";
  }

  if (
    stepIndex < currentIndex
  ) {
    return "complete";
  }

  if (
    stepIndex === currentIndex
  ) {
    return "active";
  }

  return "pending";
}


function StepNode({
  state,
  icon: Icon,
}: {
  state: VisualState;
  icon: LucideIcon;
}) {
  if (state === "complete") {
    return (
      <div className="relative z-10 flex h-[46px] w-[46px] items-center justify-center rounded-full border border-emerald-400/25 bg-[#0b2826] shadow-[0_0_24px_rgba(52,211,153,0.10)]">
        <Check
          size={20}
          strokeWidth={2.4}
          className="text-emerald-400"
        />
      </div>
    );
  }

  if (state === "active") {
    return (
      <div className="relative z-10">
        <div className="absolute inset-[-6px] animate-pulse rounded-full border border-cyan-400/20 bg-cyan-400/[0.04]" />

        <div className="relative flex h-[46px] w-[46px] items-center justify-center rounded-full border border-cyan-400/35 bg-[#0a2530] shadow-[0_0_28px_rgba(34,211,238,0.16)]">
          <LoaderCircle
            size={21}
            strokeWidth={2}
            className="animate-spin text-cyan-300"
          />
        </div>
      </div>
    );
  }

  return (
    <div className="relative z-10 flex h-[46px] w-[46px] items-center justify-center rounded-full border border-white/[0.08] bg-[#0b1927]">
      <Icon
        size={19}
        strokeWidth={1.7}
        className="text-slate-600"
      />
    </div>
  );
}


export default function RunPipeline({
  run,
}: RunPipelineProps) {
  const active =
    run.status === "running";

  const now =
    useCurrentTime(
      active,
    );

  const currentIndex =
    PIPELINE_STEPS.findIndex(
      (step) =>
        step.stage === run.stage,
    );

  const progress = useMemo(
    () => {
      if (
        run.status === "completed"
      ) {
        return 100;
      }

      if (
        run.status !== "running" ||
        currentIndex < 0
      ) {
        return 0;
      }

      if (
        PIPELINE_STEPS.length <= 1
      ) {
        return 0;
      }

      return (
        currentIndex /
        (
          PIPELINE_STEPS.length -
          1
        )
      ) * 100;
    },
    [
      run.status,
      currentIndex,
    ],
  );

  const duration =
    formatDuration(
      run.started_at_utc,
      run.finished_at_utc,
      now,
    );

  return (
    <section className="relative overflow-hidden rounded-[16px] border border-white/[0.07] bg-[#091725]/95 shadow-[0_18px_55px_rgba(0,0,0,0.16)]">
      <div className="pointer-events-none absolute inset-x-0 top-0 h-px bg-gradient-to-r from-transparent via-cyan-300/[0.12] to-transparent" />

      <div className="flex min-h-[70px] items-center justify-between border-b border-white/[0.06] px-5">
        <div>
          <div className="flex items-center gap-2.5">
            <h2 className="text-[16px] font-semibold tracking-[-0.02em] text-slate-100">
              Run Pipeline
            </h2>

            {run.status ===
              "running" && (
              <div className="flex items-center gap-1.5 rounded-full border border-cyan-400/10 bg-cyan-400/[0.04] px-2 py-1">
                <span className="h-1.5 w-1.5 animate-pulse rounded-full bg-cyan-400 shadow-[0_0_8px_rgba(34,211,238,0.7)]" />

                <span className="text-[10px] font-semibold uppercase tracking-[0.08em] text-cyan-300">
                  Live
                </span>
              </div>
            )}
          </div>

          <p className="mt-1 text-[12px] text-slate-500">
            Scientific workflow execution
            and independent verification
          </p>
        </div>

        <div className="flex items-center gap-3">
          <div
            className={[
              "flex h-8 items-center gap-2 rounded-lg border px-3 text-[12px] font-semibold",
              statusClasses(
                run,
              ),
            ].join(" ")}
          >
            {run.status ===
            "failed" ? (
              <XCircle
                size={14}
              />
            ) : (
              <Circle
                size={8}
                fill="currentColor"
              />
            )}

            {statusLabel(
              run,
            )}
          </div>

          <div className="flex h-8 items-center gap-2 rounded-lg border border-white/[0.06] bg-white/[0.02] px-3">
            <Clock3
              size={14}
              className="text-slate-500"
            />

            <span className="font-mono text-[13px] font-medium tabular-nums text-slate-300">
              {duration}
            </span>
          </div>
        </div>
      </div>

      <div className="px-6 pb-6 pt-7">
        <div className="relative">
          <div className="absolute left-[8.3%] right-[8.3%] top-[23px] hidden h-px bg-white/[0.07] lg:block">
            <div
              className="h-full bg-gradient-to-r from-emerald-400 via-cyan-400 to-cyan-300 transition-[width] duration-700 ease-out"
              style={{
                width:
                  `${progress}%`,
              }}
            />
          </div>

          <div className="grid grid-cols-2 gap-x-3 gap-y-6 md:grid-cols-3 lg:grid-cols-6">
            {PIPELINE_STEPS.map(
              (
                step,
                index,
              ) => {
                const visualState =
                  getVisualState(
                    run,
                    index,
                    currentIndex,
                  );

                const Icon =
                  step.icon;

                return (
                  <div
                    key={
                      step.stage
                    }
                    className="relative flex min-w-0 flex-col items-center text-center"
                  >
                    <StepNode
                      state={
                        visualState
                      }
                      icon={Icon}
                    />

                    <p
                      className={[
                        "mt-3 text-[13px] font-semibold tracking-[-0.01em]",
                        visualState ===
                        "complete"
                          ? "text-slate-200"
                          : visualState ===
                              "active"
                            ? "text-cyan-200"
                            : "text-slate-500",
                      ].join(" ")}
                    >
                      {step.label}
                    </p>

                    <p className="mt-1 text-[11px] text-slate-600">
                      {
                        step.description
                      }
                    </p>

                    <div className="mt-2">
                      {visualState ===
                        "complete" && (
                        <span className="text-[10px] font-medium uppercase tracking-[0.08em] text-emerald-400/80">
                          Complete
                        </span>
                      )}

                      {visualState ===
                        "active" && (
                        <span className="text-[10px] font-medium uppercase tracking-[0.08em] text-cyan-300">
                          Running
                        </span>
                      )}

                      {visualState ===
                        "pending" && (
                        <span className="text-[10px] font-medium uppercase tracking-[0.08em] text-slate-700">
                          Waiting
                        </span>
                      )}
                    </div>
                  </div>
                );
              },
            )}
          </div>
        </div>

        {run.status ===
          "failed" &&
          run.error && (
            <div className="mt-6 flex items-start gap-3 rounded-xl border border-red-400/15 bg-red-400/[0.04] px-4 py-3">
              <XCircle
                size={17}
                className="mt-0.5 shrink-0 text-red-400"
              />

              <div>
                <p className="text-[12px] font-semibold text-red-300">
                  Scientific run failed
                </p>

                <p className="mt-1 text-[12px] leading-5 text-red-200/60">
                  {run.error}
                </p>
              </div>
            </div>
          )}

        <div className="mt-6 flex items-center justify-between border-t border-white/[0.05] pt-4">
          <div className="min-w-0">
            <span className="text-[11px] text-slate-600">
              Run ID
            </span>

            <span className="ml-2 font-mono text-[11px] text-slate-500">
              {run.run_id}
            </span>
          </div>

          <div className="text-[11px] text-slate-600">
            Pipeline v
            {run.pipeline_version ??
              "0.1.0"}
          </div>
        </div>
      </div>
    </section>
  );
}
