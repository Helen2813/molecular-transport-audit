import {
  Activity,
} from "lucide-react";

import {
  CartesianGrid,
  Line,
  LineChart,
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


interface RandomControlChartProps {
  modules: ModuleResult[];
}


function formatLogTick(
  value: number,
): string {
  if (value === 1) {
    return "10⁰";
  }

  if (value === 0.1) {
    return "10⁻¹";
  }

  if (value === 0.01) {
    return "10⁻²";
  }

  if (value === 0.001) {
    return "10⁻³";
  }

  if (value === 0.0001) {
    return "10⁻⁴";
  }

  return "";
}


export default function RandomControlChart({
  modules,
}: RandomControlChartProps) {
  const data = modules.map(
    (module) => ({
      module: module.module,

      empiricalP:
        module.random_control
          .empirical_p ?? 1,
    }),
  );

  return (
    <ChartCard
      title="Random Panel Specificity"
      subtitle="Empirical p · lower is better"
      icon={Activity}
      action={
        <div className="flex items-center gap-2">
          <span className="h-2 w-2 rounded-full bg-[#8b70f8] shadow-[0_0_10px_rgba(139,112,248,0.6)]" />

          <span className="text-[11px] text-slate-500">
            empirical p
          </span>
        </div>
      }
    >
      <div className="h-[285px] w-full">
        <ResponsiveContainer
          width="100%"
          height="100%"
        >
          <LineChart
            data={data}
            margin={{
              top: 12,
              right: 18,
              bottom: 8,
              left: 3,
            }}
          >
            <defs>
              <filter
                id="purpleGlow"
                x="-50%"
                y="-50%"
                width="200%"
                height="200%"
              >
                <feGaussianBlur
                  stdDeviation="3"
                  result="blur"
                />

                <feMerge>
                  <feMergeNode
                    in="blur"
                  />

                  <feMergeNode
                    in="SourceGraphic"
                  />
                </feMerge>
              </filter>
            </defs>

            <CartesianGrid
              vertical={false}
              stroke="rgba(148,163,184,0.08)"
              strokeDasharray="4 6"
            />

            <XAxis
              dataKey="module"
              tick={{
                fill: "#cbd5e1",
                fontSize: 13,
                fontWeight: 600,
              }}
              axisLine={false}
              tickLine={false}
              dy={8}
            />

            <YAxis
              scale="log"
              domain={[
                0.0001,
                1,
              ]}
              ticks={[
                1,
                0.1,
                0.01,
                0.001,
                0.0001,
              ]}
              tickFormatter={
                formatLogTick
              }
              allowDataOverflow
              tick={{
                fill: "#64748b",
                fontSize: 12,
              }}
              axisLine={false}
              tickLine={false}
              width={42}
            />

            <ReferenceLine
              y={0.01}
              stroke="#8b70f8"
              strokeOpacity={0.48}
              strokeDasharray="5 5"
              label={{
                value: "α = 0.01",
                position: "insideTopRight",
                fill: "#8b70f8",
                fontSize: 11,
              }}
            />

            <Tooltip
              content={
                <ChartTooltip
                  formatter={(
                    value,
                  ) =>
                    value.toFixed(6)
                  }
                />
              }
            />

            <Line
              type="monotone"
              dataKey="empiricalP"
              name="Empirical p"
              stroke="#8b70f8"
              strokeWidth={2.5}
              activeDot={{
                r: 6,
                fill: "#a78bfa",
                stroke: "#ddd6fe",
                strokeWidth: 2,
              }}
              dot={{
                r: 4.5,
                fill: "#8b70f8",
                stroke: "#c4b5fd",
                strokeWidth: 1.5,
                filter:
                  "url(#purpleGlow)",
              }}
              isAnimationActive
            />
          </LineChart>
        </ResponsiveContainer>
      </div>
    </ChartCard>
  );
}
