import {
  CalendarDays,
  ChevronDown,
  Filter,
  Plus,
  Search,
} from "lucide-react";

export default function TopBar() {
  return (
    <header className="fixed left-[248px] right-0 top-0 z-20 h-[72px] border-b border-white/[0.07] bg-[#07131f]/90 backdrop-blur-xl">
      <div className="flex h-full items-center justify-between px-6">
        <div className="flex h-11 w-[420px] items-center gap-3 rounded-xl border border-white/[0.08] bg-[#0c1a28] px-4 shadow-[inset_0_1px_0_rgba(255,255,255,0.02)]">
          <Search
            size={18}
            className="shrink-0 text-slate-500"
          />

          <input
            type="text"
            placeholder="Search runs, modules, cohorts..."
            className="min-w-0 flex-1 bg-transparent text-[13px] text-slate-200 outline-none placeholder:text-slate-500"
          />

          <div className="rounded-md border border-white/[0.06] bg-white/[0.025] px-1.5 py-0.5 text-[11px] text-slate-500">
            ⌘K
          </div>
        </div>

        <div className="flex items-center gap-3">
          <button
            type="button"
            className="flex h-11 items-center gap-3 rounded-xl border border-white/[0.09] bg-[#091725] px-4 text-[13px] font-medium text-slate-200 transition hover:border-white/[0.15] hover:bg-[#0d1c2b]"
          >
            <CalendarDays
              size={18}
              className="text-slate-300"
            />

            May 12 – May 19, 2025

            <ChevronDown
              size={15}
              className="ml-1 text-slate-500"
            />
          </button>

          <button
            type="button"
            className="flex h-11 items-center gap-3 rounded-xl border border-white/[0.09] bg-[#091725] px-4 text-[13px] font-medium text-slate-200 transition hover:border-white/[0.15] hover:bg-[#0d1c2b]"
          >
            <Filter
              size={18}
              className="text-slate-300"
            />

            All Statuses

            <ChevronDown
              size={15}
              className="ml-1 text-slate-500"
            />
          </button>

          <button
            type="button"
            className="flex h-11 items-center gap-2.5 rounded-xl border border-blue-400/20 bg-gradient-to-r from-[#245cf5] to-[#7146ef] px-5 text-[13px] font-semibold text-white shadow-[0_8px_30px_rgba(50,80,255,0.18)] transition hover:brightness-110"
          >
            <Plus
              size={18}
              strokeWidth={2.2}
            />

            New Run
          </button>

          <button
            type="button"
            className="ml-1 flex h-11 items-center gap-3 rounded-xl border border-white/[0.08] bg-[#091725] px-2.5"
            aria-label="Open user menu"
          >
            <div className="relative flex h-8 w-8 items-center justify-center rounded-full bg-gradient-to-br from-slate-500 to-slate-700 text-xs font-semibold text-white">
              OA

              <span className="absolute -bottom-0.5 -right-0.5 h-2.5 w-2.5 rounded-full border-2 border-[#091725] bg-emerald-400" />
            </div>

            <ChevronDown
              size={15}
              className="text-slate-500"
            />
          </button>
        </div>
      </div>
    </header>
  );
}
