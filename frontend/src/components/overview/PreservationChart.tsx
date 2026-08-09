import {
  BarChart3,
} from "lucide-react";

import {
  Bar,
  BarChart,
  CartesianGrid,
  LabelList,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

import type {
  ModuleResult,
} from "../../types/audit";

import ChartCard from "./ChartCard";
import ChartTooltip from "./ChartTooltip";


interface PreservationChartProps {
  modules: ModuleResult[];
}


export default function PreservationChart({
  modules,
}: PreservationChartProps) {
  const data = modules.map(
    (module) => ({
      module: module.module,

      edge:
        module.direct_preservation
          .edge_spearman ?? 0,

      loading:
        module.direct_preservation
          .loading_spearman ?? 0,
    }),
  );

  return (
    <ChartCard
      title="Preservation by Module"
      subtitle="Edge and loading concordance"
      icon={BarChart3}
      action={
        <div className="hidden items-center gap-4 xl:flex">
          <div className="flex items-center gap-2">
            <span className="h-2 w-2 rounded-full bg-[#25c8d9]" />

            <span className="text-[11px] text-slate-500">
              Edge
            </span>
          </div>

          <div className="flex items-center gap-2">
            <span className="h-2 w-2 rounded-full bg-[#8068f5]" />

            <span className="text-[11px] text-slate-500">
              Loading
            </span>
          </div>
        </div>
      }
    >
      <div className="h-[285px] w-full">
        <ResponsiveContainer
          width="100%"
          height="100%"
        >
          <BarChart
            data={data}
            layout="vertical"
            margin={{
              top: 4,
              right: 42,
              bottom: 8,
              left: 2,
            }}
            barCategoryGap={12}
            barGap={3}
          >
            <CartesianGrid
              horizontal={false}
              stroke="rgba(148,163,184,0.08)"
              strokeDasharray="4 6"
            />

            <XAxis
              type="number"
              domain={[-1, 1]}
              ticks={[
                -1,
                -0.5,
                0,
                0.5,
                1,
              ]}
              tick={{
                fill: "#64748b",
                fontSize: 12,
              }}
              axisLine={false}
              tickLine={false}
            />

            <YAxis
              type="category"
              dataKey="module"
              width={42}
              tick={{
                fill: "#cbd5e1",
                fontSize: 13,
                fontWeight: 600,
              }}
              axisLine={false}
              tickLine={false}
            />

            <ReferenceLine
              x={0}
              stroke="rgba(148,163,184,0.22)"
            />

            <Tooltip
              cursor={{
                fill:
                  "rgba(255,255,255,0.025)",
              }}
              content={
                <ChartTooltip />
              }
            />

            <Bar
              dataKey="edge"
              name="Edge Spearman ρ"
              fill="#25c8d9"
              radius={[0, 5, 5, 0]}
              maxBarSize={15}
              isAnimationActive
            >
              <LabelList
                dataKey="edge"
                position="right"
                formatter={(
                  value: number,
                ) =>
                  value.toFixed(3)
                }
                fill="#94a3b8"
                fontSize={11}
              />
            </Bar>

            <Bar
              dataKey="loading"
              name="Loading Spearman ρ"
              fill="#8068f5"
              radius={[0, 5, 5, 0]}
              maxBarSize={15}
              isAnimationActive
            />
          </BarChart>
        </ResponsiveContainer>
      </div>
    </ChartCard>
  );
}
