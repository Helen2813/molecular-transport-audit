import type {
  ModuleResult,
} from "../../types/audit";

import PreservationChart from "./PreservationChart";
import RandomControlChart from "./RandomControlChart";
import ReliabilityChart from "./ReliabilityChart";


interface OverviewChartsProps {
  modules: ModuleResult[];
}


export default function OverviewCharts({
  modules,
}: OverviewChartsProps) {
  return (
    <section className="grid grid-cols-1 gap-4 xl:grid-cols-2 2xl:grid-cols-[1.15fr_1fr_1fr]">
      <PreservationChart
        modules={modules}
      />

      <ReliabilityChart
        modules={modules}
      />

      <div className="xl:col-span-2 2xl:col-span-1">
        <RandomControlChart
          modules={modules}
        />
      </div>
    </section>
  );
}
