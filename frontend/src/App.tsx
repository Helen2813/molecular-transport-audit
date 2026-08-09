import Sidebar from "./components/layout/Sidebar";
import TopBar from "./components/layout/TopBar";
import OverviewPage from "./pages/OverviewPage";

export default function App() {
  return (
    <div className="min-h-screen bg-transparent text-white">
      <Sidebar />

      <TopBar />

      <OverviewPage />
    </div>
  );
}
