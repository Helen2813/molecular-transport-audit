interface TooltipEntry {
  name?: string;
  value?: number | string;
  color?: string;
}


interface ChartTooltipProps {
  active?: boolean;
  label?: string;
  payload?: TooltipEntry[];
  formatter?: (
    value: number,
  ) => string;
}


export default function ChartTooltip({
  active,
  label,
  payload,
  formatter = (value) =>
    value.toFixed(3),
}: ChartTooltipProps) {
  if (
    !active ||
    !payload ||
    payload.length === 0
  ) {
    return null;
  }

  return (
    <div className="min-w-[160px] rounded-xl border border-white/[0.10] bg-[#07131f]/95 px-3.5 py-3 shadow-[0_16px_50px_rgba(0,0,0,0.45)] backdrop-blur-xl">
      <p className="mb-2.5 text-[13px] font-semibold text-white">
        {label}
      </p>

      <div className="space-y-2">
        {payload.map(
          (entry, index) => {
            const numeric =
              typeof entry.value ===
              "number"
                ? entry.value
                : Number(entry.value);

            return (
              <div
                key={`${entry.name}-${index}`}
                className="flex items-center justify-between gap-5"
              >
                <div className="flex items-center gap-2">
                  <span
                    className="h-2 w-2 rounded-full"
                    style={{
                      backgroundColor:
                        entry.color,
                    }}
                  />

                  <span className="text-[12px] text-slate-400">
                    {entry.name}
                  </span>
                </div>

                <span className="font-mono text-[12px] font-medium tabular-nums text-slate-100">
                  {Number.isFinite(
                    numeric,
                  )
                    ? formatter(
                        numeric,
                      )
                    : "—"}
                </span>
              </div>
            );
          },
        )}
      </div>
    </div>
  );
}
