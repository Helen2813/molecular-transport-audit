export default function OverviewSkeleton() {
  return (
    <div className="animate-pulse">
      <div className="h-8 w-[360px] rounded-lg bg-white/[0.05]" />

      <div className="mt-3 h-4 w-[210px] rounded bg-white/[0.035]" />

      <div className="mt-6 grid grid-cols-1 gap-3 sm:grid-cols-2 xl:grid-cols-5">
        {Array.from({
          length: 5,
        }).map((_, index) => (
          <div
            key={index}
            className="h-[100px] rounded-[13px] border border-white/[0.05] bg-white/[0.025]"
          />
        ))}
      </div>
    </div>
  );
}
