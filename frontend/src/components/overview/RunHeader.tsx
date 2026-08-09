import {
  ExternalLink,
  Star,
} from "lucide-react";

interface RunHeaderProps {
  externalCohort: string;
}

export default function RunHeader({
  externalCohort,
}: RunHeaderProps) {
  return (
    <div className="mb-5 flex items-start justify-between gap-6">
      <div>
        <div className="flex items-center gap-2.5">
          <h1 className="text-[25px] font-semibold tracking-[-0.03em] text-[#f5f7fb]">
            {externalCohort} External Canine Audit
          </h1>

          <button
            type="button"
            className="rounded-md p-1 text-slate-500 transition hover:bg-white/[0.04] hover:text-slate-300"
            aria-label="Favorite run"
          >
            <Star
              size={17}
              strokeWidth={1.7}
            />
          </button>
        </div>

        <p className="mt-1.5 text-[13px] text-slate-400">
          Validated transport audit overview
        </p>
      </div>

      <button
        type="button"
        className="flex h-10 shrink-0 items-center gap-2 rounded-lg border border-white/[0.08] bg-[#0a1825] px-4 text-[12px] font-medium text-slate-200 transition hover:border-white/[0.14] hover:bg-[#0d1d2c]"
      >
        View Run Details

        <ExternalLink
          size={14}
          className="text-slate-500"
        />
      </button>
    </div>
  );
}
