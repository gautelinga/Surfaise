import dolfin as df
import surfaise as sf
import os
import ufl
from surfaise.common.io import (
    Timeseries, save_checkpoint, load_checkpoint,
    load_parameters)
from surfaise.common.cmd import (
    mpi_max, parse_command_line, info_blue, info_cyan)


parameters = dict(
    R=30.,  # Major Radius. "Default" values: R=30, r=10
    r=10.,  # Minir Radius.
    res=100,  # Resolution
    dt=1e-1,
    rho=0.01,
    mu=1.0,
    restart_folder=None,
    folder="results_ns_torus",
    t_0=0.0,
    tstep=0,
    T=2000,
    dump_intv=5,
    stats_intv=1,
    checkpoint_intv=50,
    verbose=True,
    init_mode="random",
    #init_mode="nothing",
    alpha=0.0,
)
cmd_kwargs = sf.cmd.parse_command_line()
parameters.update(**cmd_kwargs)
if parameters["restart_folder"]:
    sf.io.load_parameters(parameters, os.path.join(
        parameters["restart_folder"], "parameters.dat"))
    parameters.update(**cmd_kwargs)

R = parameters["R"]
r = parameters["r"]
res = parameters["res"]
dt = sf.TimeStepSelector(parameters["dt"])
rho = df.Constant(parameters["rho"])
mu = df.Constant(parameters["mu"])

geo_map = sf.TorusMap(R, r)
#geo_map = sf.EllipsoidMap(R, R, R, verbose=True, eps=1e-1)
geo_map.initialize(res, restart_folder=parameters["restart_folder"])

W = geo_map.mixed_space((geo_map.ref_vel, geo_map.ref_el))

# W = geo_map.mixed_space((geo_map.ref_el, geo_map.ref_el, geo_map.ref_el))
# field_names = ("psi", "mu", "nu")

# Define trial and test functions
du = df.TrialFunction(W)
v, q = df.TestFunctions(W)

# Define functions
w = df.TrialFunction(W)
w_ = df.Function(W, name="u_")  # current solution
w_1 = df.Function(W, name="u_1")  # solution from previous converged step

# Split mixed functions
u,  p = df.split(w)
u_, p_ = df.split(w_)
u_1, p_1 = df.split(w_1)

# Create intial conditions
if parameters["restart_folder"] is None:
    init_mode = parameters["init_mode"]
    if init_mode == "random":
        w_init = sf.ics.RandomIC(w_, degree=1)
    else:
        exit("Unknown IC")
    w_1.interpolate(w_init)
    w_.assign(w_1)
else:
    sf.io.load_checkpoint(parameters["restart_folder"], w_, w_1)


f = df.Expression(("0.1*exp(-pow(x[1]-1.57,2)/2*0.01)", "0."), degree=2)

# Define some UFL indices:
i, j, k, l = ufl.Index(), ufl.Index(), ufl.Index(), ufl.Index()

m_NS = (rho / dt * geo_map.g_ab[i, j] * (u_[i]-u_1[i]) * v[j]
        + rho * geo_map.g_ab[i, j] * u_[k] *
        geo_map.CovD10(u_)[k, i] * v[j]
        + mu * geo_map.g_ab[i, k] * geo_map.gab[j, l] *
        geo_map.CovD10(u_)[i, j] * geo_map.CovD10(v)[k, l]
        + mu * geo_map.K * geo_map.g_ab[i, j] * u_[i] * v[j]
        - p_ * geo_map.CovD10(v)[i, i]
        + q * geo_map.CovD10(u_)[i, i]
        - f[i]*v[i])
F = geo_map.form(m_NS)

J = df.derivative(F, w_, du=w)

#problem = df.NonlinearVariationalProblem(F, w_, J=J, bcs=[u_bc])
problem = df.NonlinearVariationalProblem(F, w_, J=J)
solver = df.NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-8
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-5
solver.parameters["newton_solver"]["maximum_iterations"] = 16
#solver.parameters["newton_solver"]["linear_solver"] = "gmres"
#solver.parameters["newton_solver"]["preconditioner"] = "default"
# solver.parameters["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
# solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-8
# solver.parameters["newton_solver"]["krylov_solver"]["monitor_convergence"] = False
#solver.parameters["newton_solver"]["krylov_solver"]["maximum_iterations"] = 1000

# solver.parameters["linear_solver"] = "gmres"
# solver.parameters["preconditioner"] = "jacobi"

df.parameters["form_compiler"]["optimize"] = True
df.parameters["form_compiler"]["cpp_optimize"] = True

#
t = parameters["t_0"]
tstep = parameters["tstep"]
T = parameters["T"]

# Output file
ts = Timeseries(parameters["folder"], w_,
                ("u", "p"), geo_map, tstep,
                parameters=parameters,
                restart_folder=parameters["restart_folder"])

E_kin = 0.5*rho*geo_map.g_ab[i, j]*u_[i]*u_[j]
divu = geo_map.CovD10(u_)[i, i]
U = [sum([geo_map.get_function(xi + "_," + vj)*u_[dj]
          for dj, vj in enumerate(geo_map.AXIS_REF)])
     for xi in geo_map.AXIS]

ts.add_field(E_kin, "E_kin")
ts.add_field(divu, "divu")
ts.add_field(U, "U")

# Step in time
ts.dump(tstep)

while t < T:
    tstep += 1
    info_cyan("tstep = {}, time = {}".format(tstep, t))

    w_1.assign(w_)

    converged = False
    while not converged:
        try:
            solver.solve()
            converged = True
        except:
            info_blue("Did not converge. Chopping timestep.")
            dt.chop()
            info_blue("New timestep is: dt = {}".format(dt.get()))

    # Update time with final dt value
    t += dt.get()

    if tstep % 1 == 0:
        ts.dump(t)
        # Assigning timestep size according to grad_mu_max:
        # dt_prev = dt.get()
        # dt.set(min(0.05/grad_mu_max, T-t))
        info_blue("dt = {}".format(dt.get()))
        # ts.dump_stats(t, "data")

    if tstep % parameters["checkpoint_intv"] == 0 or t >= T:
        save_checkpoint(tstep, t, geo_map.ref_mesh,
                        w_, w_1, ts.folder, parameters)
