import {
  ShieldCheck,
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


interface ReliabilityChartProps {
  modules: ModuleResult[];
}


export default function ReliabilityChart({
  modules,
}: ReliabilityChartProps) {
  const data = modules.map(
    (module) => ({
      module: module.module,

      reliability:
        module.reliability
          .split_half_median ?? 0,
    }),
  );

  return (
    <ChartCard
      title="Reliability"
      subtitle="Split-half median ρ"
      icon={ShieldCheck}
      action={
        <div className="rounded-md border border-cyan-400/10 bg-cyan-400/[0.05] px-2 py-1 text-[10px] font-medium text-cyan-300">
          threshold 0.70
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
              right: 44,
              bottom: 8,
              left: 2,
            }}
            barCategoryGap={18}
          >
            <CartesianGrid
              horizontal={false}
              stroke="rgba(148,163,184,0.08)"
              strokeDasharray="4 6"
            />

            <XAxis
              type="number"
              domain={[0, 1]}
              ticks={[
                0,
                0.25,
                0.5,
                0.75,
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
              x={0.7}
              stroke="#22d3ee"
              strokeOpacity={0.45}
              strokeDasharray="5 5"
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
              dataKey="reliability"
              name="Split-half median ρ"
              fill="#22b8cf"
              radius={[0, 6, 6, 0]}
              maxBarSize={22}
              isAnimationActive
            >
              <LabelList
                dataKey="reliability"
                position="right"
                formatter={(
                  value: number,
                ) =>
                  value.toFixed(3)
                }
                fill="#cbd5e1"
                fontSize={12}
                fontWeight={500}
              />
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </div>
    </ChartCard>
  );
}
