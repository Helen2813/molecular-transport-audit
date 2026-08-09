import {
  Activity,
  Box,
  Braces,
  ChevronLeft,
  CirclePlay,
  FileText,
  Gauge,
  House,
  Settings,
  ShieldCheck,
  Shuffle,
} from "lucide-react";

const navigation = [
  {
    label: "Overview",
    icon: House,
    active: true,
  },
  {
    label: "Runs",
    icon: CirclePlay,
  },
  {
    label: "Modules",
    icon: Box,
  },
  {
    label: "Reliability",
    icon: ShieldCheck,
  },
  {
    label: "Random Controls",
    icon: Shuffle,
  },
  {
    label: "Validation",
    icon: ShieldCheck,
  },
  {
    label: "Reports",
    icon: FileText,
  },
  {
    label: "API",
    icon: Braces,
  },
  {
    label: "Settings",
    icon: Settings,
  },
];

function LogoMark() {
  return (
    <div className="relative h-12 w-12">
      <div className="absolute left-1/2 top-1/2 h-8 w-8 -translate-x-1/2 -translate-y-1/2 rotate-45 rounded-[8px] border border-cyan-400/70" />

      <span className="absolute left-[5px] top-[7px] h-2.5 w-2.5 rounded-full bg-cyan-400 shadow-[0_0_14px_rgba(34,211,238,0.55)]" />

      <span className="absolute right-[5px] top-[15px] h-2.5 w-2.5 rounded-full bg-blue-500" />

      <span className="absolute bottom-[6px] left-[17px] h-2.5 w-2.5 rounded-full bg-cyan-500" />
    </div>
  );
}

export default function Sidebar() {
  return (
    <aside className="fixed inset-y-0 left-0 z-30 flex w-[248px] flex-col border-r border-white/[0.07] bg-[#06111d]/95">
      <div className="flex h-[92px] items-center gap-2.5 px-6">
        <LogoMark />

        <div className="leading-[1.14]">
          <div className="text-[17px] font-semibold tracking-[-0.02em] text-white">
            Molecular
          </div>

          <div className="text-[17px] font-semibold tracking-[-0.02em] text-white">
            Transport Audit
          </div>
        </div>
      </div>

      <nav className="mt-2 flex flex-1 flex-col gap-1.5 px-4">
        {navigation.map((item) => {
          const Icon = item.icon;

          return (
            <button
              key={item.label}
              type="button"
              className={[
                "group flex h-12 w-full items-center gap-4 rounded-xl px-4 text-left text-[14px] font-medium transition",
                item.active
                  ? "border border-cyan-400/20 bg-gradient-to-r from-blue-500/15 to-cyan-400/[0.07] text-[#41a6ff] shadow-[inset_0_0_24px_rgba(24,119,242,0.05)]"
                  : "border border-transparent text-slate-300 hover:bg-white/[0.035] hover:text-white",
              ].join(" ")}
            >
              <Icon
                size={20}
                strokeWidth={1.8}
                className={
                  item.active
                    ? "text-[#41a6ff]"
                    : "text-slate-400 group-hover:text-slate-200"
                }
              />

              {item.label}
            </button>
          );
        })}
      </nav>

      <div className="px-4 pb-4">
        <div className="mb-5 rounded-xl border border-white/[0.08] bg-white/[0.015] px-4 py-3.5">
          <div className="flex items-center gap-2 text-xs text-slate-300">
            <span className="h-2 w-2 rounded-full bg-emerald-400 shadow-[0_0_10px_rgba(52,211,153,0.55)]" />

            System Status
          </div>

          <div className="mt-1.5 text-xs text-emerald-400">
            All systems operational
          </div>
        </div>

        <button
          type="button"
          className="flex h-10 w-full items-center gap-3 px-3 text-sm text-slate-400 transition hover:text-white"
        >
          <ChevronLeft size={17} />

          Collapse
        </button>
      </div>
    </aside>
  );
}
