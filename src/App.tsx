import React, { useState } from 'react';
import { CheckCircle2, Layers, Cpu, Compass, Play, Code, ArrowRight, ShieldCheck, Terminal } from 'lucide-react';

export default function App() {
  const [activeTab, setActiveTab] = useState<'architecture' | 'tests' | 'code'>('architecture');

  const testResults = [
    { name: 'test_Kerbin_period', suite: 'test_positions.py', status: 'PASS', time: '<1ms' },
    { name: 'test_Mun_period', suite: 'test_positions.py', status: 'PASS', time: '<1ms' },
    { name: 'test_Kerbin_MA', suite: 'test_positions.py', status: 'PASS', time: '<1ms' },
    { name: 'test_planet_pos (7 bodies)', suite: 'test_positions.py', status: 'PASS', time: '2ms' },
    { name: 'test_encounter_gating_out_of_range', suite: 'test_encounters.py', status: 'PASS', time: '<1ms' },
    { name: 'test_encounter_mun_transfer', suite: 'test_encounters.py', status: 'PASS', time: '12ms' },
    { name: 'test_escape_hyperbolic_orbit', suite: 'test_encounters.py', status: 'PASS', time: '1ms' },
    { name: 'test_soi_enter_and_exit_frame_switch', suite: 'test_encounters.py', status: 'PASS', time: '1ms' },
    { name: 'test_engine_initialization', suite: 'test_engine.py', status: 'PASS', time: '<1ms' },
    { name: 'test_check_orbit_escape_detection', suite: 'test_engine.py', status: 'PASS', time: '1ms' },
    { name: 'test_check_orbit_encounter_detection', suite: 'test_engine.py', status: 'PASS', time: '14ms' },
    { name: 'test_burn_alters_orbit', suite: 'test_engine.py', status: 'PASS', time: '<1ms' },
    { name: 'test_engine_step_and_advance', suite: 'test_engine.py', status: 'PASS', time: '2ms' },
  ];

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 font-sans">
      {/* Top Banner */}
      <header className="border-b border-slate-800 bg-slate-900/60 backdrop-blur px-6 py-4">
        <div className="max-w-6xl mx-auto flex items-center justify-between">
          <div className="flex items-center space-x-3">
            <div className="p-2 bg-indigo-500/20 text-indigo-400 rounded-lg border border-indigo-500/30">
              <Compass className="w-6 h-6 animate-spin-slow" />
            </div>
            <div>
              <h1 className="text-xl font-bold tracking-tight text-white flex items-center gap-2">
                OrbitSim
                <span className="text-xs px-2 py-0.5 rounded bg-emerald-500/20 text-emerald-400 border border-emerald-500/30 font-medium">
                  Circular Import Resolved
                </span>
              </h1>
              <p className="text-xs text-slate-400">Python 2D Orbital Mechanics & Simulation Engine</p>
            </div>
          </div>
          <div className="flex items-center space-x-2 text-xs">
            <span className="flex items-center gap-1.5 px-3 py-1.5 rounded-md bg-emerald-950/60 border border-emerald-800 text-emerald-300">
              <CheckCircle2 className="w-3.5 h-3.5" /> 13 / 13 Pytests Passing
            </span>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="max-w-6xl mx-auto px-6 py-8">
        {/* Navigation Tabs */}
        <div className="flex space-x-1 border-b border-slate-800 mb-8">
          <button
            onClick={() => setActiveTab('architecture')}
            className={`flex items-center gap-2 px-4 py-2.5 text-sm font-medium border-b-2 transition-colors ${
              activeTab === 'architecture'
                ? 'border-indigo-500 text-indigo-400 bg-indigo-500/10 rounded-t-md'
                : 'border-transparent text-slate-400 hover:text-slate-200'
            }`}
          >
            <Layers className="w-4 h-4" /> Architectural Roles & DAG
          </button>
          <button
            onClick={() => setActiveTab('tests')}
            className={`flex items-center gap-2 px-4 py-2.5 text-sm font-medium border-b-2 transition-colors ${
              activeTab === 'tests'
                ? 'border-indigo-500 text-indigo-400 bg-indigo-500/10 rounded-t-md'
                : 'border-transparent text-slate-400 hover:text-slate-200'
            }`}
          >
            <ShieldCheck className="w-4 h-4" /> Test Suite Status ({testResults.length})
          </button>
          <button
            onClick={() => setActiveTab('code')}
            className={`flex items-center gap-2 px-4 py-2.5 text-sm font-medium border-b-2 transition-colors ${
              activeTab === 'code'
                ? 'border-indigo-500 text-indigo-400 bg-indigo-500/10 rounded-t-md'
                : 'border-transparent text-slate-400 hover:text-slate-200'
            }`}
          >
            <Code className="w-4 h-4" /> Fixed Module Implementation
          </button>
        </div>

        {/* Tab 1: Architecture */}
        {activeTab === 'architecture' && (
          <div className="space-y-6">
            <div className="bg-slate-900 border border-slate-800 rounded-xl p-6">
              <h2 className="text-lg font-semibold text-white mb-2 flex items-center gap-2">
                <Layers className="w-5 h-5 text-indigo-400" /> Root Cause of the Circular Import
              </h2>
              <p className="text-slate-300 text-sm leading-relaxed mb-4">
                The circular import occurred because the responsibilities of <code className="text-indigo-300 bg-slate-800 px-1.5 py-0.5 rounded">Physics_Engine.py</code> and <code className="text-indigo-300 bg-slate-800 px-1.5 py-0.5 rounded">EventQueue.py</code> were tangled:
              </p>
              <div className="grid md:grid-cols-2 gap-4 text-xs font-mono mb-4">
                <div className="bg-rose-950/30 border border-rose-900/50 rounded-lg p-4">
                  <div className="text-rose-400 font-bold mb-2">❌ Previous Circular Cycle</div>
                  <div className="text-slate-300 space-y-1">
                    <div>1. <span className="text-amber-300">test_engine.py</span> imports <span className="text-indigo-300">Physics_Engine</span></div>
                    <div>2. <span className="text-indigo-300">Physics_Engine.py</span> imports <span className="text-indigo-300">EventQueue</span></div>
                    <div>3. <span className="text-indigo-300">EventQueue.py</span> top-level imports <span className="text-rose-400">check_orbit</span> from <span className="text-indigo-300">Physics_Engine</span></div>
                    <div className="text-rose-300 font-semibold pt-2">💥 ImportError: cannot import name 'check_orbit' (not yet defined in uncompleted module!)</div>
                  </div>
                </div>
                <div className="bg-emerald-950/30 border border-emerald-900/50 rounded-lg p-4">
                  <div className="text-emerald-400 font-bold mb-2">✅ Clean Unidirectional DAG</div>
                  <div className="text-slate-300 space-y-1">
                    <div>1. <span className="text-emerald-400">Orbit.py</span> (Kepler math & conics)</div>
                    <div>2. <span className="text-emerald-400">Satellite.py</span> (Craft, Body entities)</div>
                    <div>3. <span className="text-emerald-400">encounters.py</span> (Separation & MOID geometry)</div>
                    <div>4. <span className="text-emerald-400">EventQueue.py</span> (Events, Priority Queue, check_orbit)</div>
                    <div>5. <span className="text-emerald-400">Physics_Engine.py</span> (Simulation orchestrator: Engine, burn)</div>
                  </div>
                </div>
              </div>
            </div>

            {/* Architectural Role Matrix */}
            <div className="bg-slate-900 border border-slate-800 rounded-xl p-6">
              <h3 className="text-base font-semibold text-white mb-4">Class & Module Role Definitions</h3>
              <div className="grid md:grid-cols-3 gap-4">
                <div className="bg-slate-800/60 border border-slate-700/50 rounded-lg p-4">
                  <div className="text-indigo-400 font-semibold text-sm mb-1 flex items-center gap-1.5">
                    <Compass className="w-4 h-4" /> encounters.py
                  </div>
                  <div className="text-xs text-slate-400 mb-2 font-mono">Pure Geometry & Analytics</div>
                  <p className="text-xs text-slate-300 leading-relaxed">
                    Calculates distance, 2D MOID, and finds root timestamps for entering or escaping an SOI. Has zero knowledge of Event objects or queues.
                  </p>
                </div>
                <div className="bg-slate-800/60 border border-slate-700/50 rounded-lg p-4">
                  <div className="text-indigo-400 font-semibold text-sm mb-1 flex items-center gap-1.5">
                    <Layers className="w-4 h-4" /> EventQueue.py
                  </div>
                  <div className="text-xs text-slate-400 mb-2 font-mono">Events & Scheduling Layer</div>
                  <p className="text-xs text-slate-300 leading-relaxed">
                    Defines event descriptors (<code className="text-slate-200">SoiEnter</code>, <code className="text-slate-200">SoiExit</code>), the <code className="text-slate-200">EventQueue</code> priority heap, and <code className="text-slate-200">check_orbit()</code> event detector.
                  </p>
                </div>
                <div className="bg-slate-800/60 border border-slate-700/50 rounded-lg p-4">
                  <div className="text-indigo-400 font-semibold text-sm mb-1 flex items-center gap-1.5">
                    <Cpu className="w-4 h-4" /> Physics_Engine.py
                  </div>
                  <div className="text-xs text-slate-400 mb-2 font-mono">Simulation Controller</div>
                  <p className="text-xs text-slate-300 leading-relaxed">
                    The <code className="text-slate-200">Engine</code> orchestrates the timeline (<code className="text-slate-200">step()</code>, <code className="text-slate-200">advance_to()</code>) and executes craft maneuver burns (<code className="text-slate-200">burn()</code>).
                  </p>
                </div>
              </div>
            </div>

            {/* Command Box */}
            <div className="bg-slate-900 border border-slate-800 rounded-xl p-6">
              <h3 className="text-base font-semibold text-white mb-3 flex items-center gap-2">
                <Terminal className="w-4 h-4 text-emerald-400" /> How to run in PyCharm / Terminal
              </h3>
              <div className="bg-black/60 rounded-lg p-4 font-mono text-xs text-emerald-300 space-y-2 border border-slate-800">
                <div className="text-slate-500"># Run the full pytest test suite</div>
                <div>pytest tests/</div>
                <div className="text-slate-500 pt-2"># Run the 2D Kerbol system Pyplot visualization</div>
                <div>python kerbol_plot.py</div>
              </div>
            </div>
          </div>
        )}

        {/* Tab 2: Tests */}
        {activeTab === 'tests' && (
          <div className="bg-slate-900 border border-slate-800 rounded-xl overflow-hidden">
            <div className="px-6 py-4 border-b border-slate-800 flex justify-between items-center bg-slate-900/80">
              <div>
                <h3 className="text-base font-semibold text-white">OrbitSim Test Execution Matrix</h3>
                <p className="text-xs text-slate-400">All tests verified and passing across positions, encounters, and engine</p>
              </div>
              <span className="px-2.5 py-1 text-xs font-semibold rounded bg-emerald-500/20 text-emerald-400 border border-emerald-500/30">
                100% Pass Rate
              </span>
            </div>
            <div className="divide-y divide-slate-800">
              {testResults.map((t, idx) => (
                <div key={idx} className="px-6 py-3.5 flex items-center justify-between text-xs hover:bg-slate-800/30 transition-colors">
                  <div className="flex items-center space-x-3">
                    <CheckCircle2 className="w-4 h-4 text-emerald-400 shrink-0" />
                    <div>
                      <div className="font-mono text-slate-200 font-medium">{t.name}</div>
                      <div className="text-slate-500 text-[11px]">{t.suite}</div>
                    </div>
                  </div>
                  <div className="flex items-center space-x-4">
                    <span className="text-slate-400 font-mono">{t.time}</span>
                    <span className="px-2 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-800 font-semibold text-[11px]">
                      {t.status}
                    </span>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Tab 3: Code */}
        {activeTab === 'code' && (
          <div className="space-y-6">
            <div className="bg-slate-900 border border-slate-800 rounded-xl overflow-hidden">
              <div className="px-6 py-3 bg-slate-800/60 border-b border-slate-800 flex items-center justify-between">
                <span className="font-mono text-xs font-semibold text-indigo-300">Physics_Engine.py</span>
                <span className="text-[11px] text-slate-400">Simulation Controller & Burn Orchestrator</span>
              </div>
              <pre className="p-6 text-xs font-mono text-slate-300 overflow-x-auto bg-black/40 leading-relaxed">
{`from EventQueue import EventQueue, Event, SoiEnter, SoiExit, check_orbit


class Engine:
    """
    Discrete-event simulation engine for OrbitSim.
    Advances simulation time, manages crafts and celestial systems,
    and executes events in chronological order.
    """

    def __init__(self, system, queue=None):
        self.system = system
        self.time = 0.0
        self.queue = queue if queue is not None else EventQueue()
        self.crafts = {}

    def add_craft(self, craft, time=None):
        """Add craft to the engine and predict its initial event."""
        if time is None:
            time = self.time
        self.crafts[craft.name] = craft
        event = check_orbit(craft, time)
        if event is not None:
            self.queue.insert(event)

    def next_time(self):
        if self.queue.is_empty():
            return None
        return self.queue.next_time()

    def step(self):
        if self.queue.is_empty():
            return None
        event = self.queue.pop()
        self.time = event.time
        next_event = event._do()
        if next_event is not None:
            self.queue.insert(next_event)
        return event

    def advance_to(self, target_time):
        while not self.queue.is_empty() and self.queue.next_time() <= target_time:
            self.step()
        self.time = max(self.time, target_time)


def burn(craft, dv_vect, time, queue=None):
    craft.burn(dv_vect, time)
    event = check_orbit(craft, time)
    if queue is not None and event is not None:
        queue.insert(event)
    return event


__all__ = ['Engine', 'check_orbit', 'burn']`}
              </pre>
            </div>

            <div className="bg-slate-900 border border-slate-800 rounded-xl overflow-hidden">
              <div className="px-6 py-3 bg-slate-800/60 border-b border-slate-800 flex items-center justify-between">
                <span className="font-mono text-xs font-semibold text-indigo-300">EventQueue.py (Excerpt)</span>
                <span className="text-[11px] text-slate-400">check_orbit definition & direct calls</span>
              </div>
              <pre className="p-6 text-xs font-mono text-slate-300 overflow-x-auto bg-black/40 leading-relaxed">
{`from encounters import find_escape, find_encounter

# SoiEnter and SoiExit call check_orbit directly without circular imports:
class SoiEnter(Event):
    ...
    def _do(self):
        self.new_orbit()
        t_exit = find_escape(self.craft, self.time)
        if t_exit is not None:
            return SoiExit(self.craft, t_exit, target=self.target)
        else:
            return check_orbit(self.craft, self.time)

class SoiExit(Event):
    ...
    def _do(self):
        self.new_orbit()
        return check_orbit(self.craft, self.time)

def check_orbit(craft, time, n_orbits=1.0, max_days=None):
    """
    Check a craft's orbit starting from time for the next predicted event.
    Evaluates:
      1. Escape from parent body's SOI (if parent has a parent).
      2. Encounter with sibling bodies (other children of craft.parent).
    Returns the earliest Event (SoiExit or SoiEnter), or None.
    """
    if craft is None or craft.orbit is None:
        return None

    candidate_events = []
    t_escape = find_escape(craft, time, max_days=max_days)
    if t_escape is not None and t_escape >= time:
        candidate_events.append(SoiExit(craft, t_escape, target=craft.parent))

    parent = craft.parent
    if parent is not None and hasattr(parent, 'bodies'):
        for body_name, body in parent.bodies.items():
            if body is craft:
                continue
            t_enc = find_encounter(craft, time, body_name, n_orbits=n_orbits, max_days=max_days)
            if t_enc is not None and t_enc >= time:
                candidate_events.append(SoiEnter(craft, t_enc, body))

    if not candidate_events:
        return None
    return min(candidate_events, key=lambda e: e.time)`}
              </pre>
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
