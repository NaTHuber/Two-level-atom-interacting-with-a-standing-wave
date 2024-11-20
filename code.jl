using DifferentialEquations, ParameterizedFunctions, Plots

# Sistema de ecuaciones diferenciales 
my_ode10 = @ode_def begin
    da₁ = -im*Ω*a₁ -im*(Ω/2)*(a₅+a₄)
    da₂ = -im*(ωᵣ + Ω)*a₂ - im*(Ω/2)*(a₃+a₆)
    da₃ = -im*(ωᵣ + Ω)*a₃ - im*(Ω/2)*(a₇+a₂)
    da₄ = -im*(4*ωᵣ + Ω)*a₄ - im*(Ω/2)*(a₁+a₈)
    da₅ = -im*(4*ωᵣ + Ω)*a₅ - im*(Ω/2)*(a₉+a₁)
    da₆ = -im*(9*ωᵣ + Ω)*a₆ - im*(Ω/2)*(a₂+a₁₀)
    da₇ = -im*(9*ωᵣ + Ω)*a₇ - im*(Ω/2)*(a₃+a₁₁)
    da₈ = -im*(16*ωᵣ + Ω)*a₈ - im*(Ω/2)*(a₄+a₁₂)
    da₉ = -im*(16*ωᵣ + Ω)*a₉ - im*(Ω/2)*(a₅+a₁₃)
    da₁₀ = -im*(25*ωᵣ + Ω)*a₁₀ - im*(Ω/2)*(a₆+a₁₄)
    da₁₁ = -im*(25*ωᵣ + Ω)*a₁₁ - im*(Ω/2)*(a₇+a₁₅)
    da₁₂ = -im*(36*ωᵣ + Ω)*a₁₂ - im*(Ω/2)*(a₈+a₁₆)
    da₁₃ = -im*(36*ωᵣ + Ω)*a₁₃ - im*(Ω/2)*(a₉+a₁₇)
    da₁₄ = -im*(49*ωᵣ + Ω)*a₁₄ - im*(Ω/2)*(a₁₀+a₁₈)
    da₁₅ = -im*(49*ωᵣ + Ω)*a₁₅ - im*(Ω/2)*(a₁₁+a₁₉)
    da₁₆ = -im*(64*ωᵣ + Ω)*a₁₆ - im*(Ω/2)*(a₁₂+a₂₀)
    da₁₇ = -im*(64*ωᵣ + Ω)*a₁₇ - im*(Ω/2)*(a₁₃+a₂₁)
    da₁₈ = -im*(81*ωᵣ + Ω)*a₁₈ - im*(Ω/2)*(a₁₄)
    da₁₉ = -im*(81*ωᵣ + Ω)*a₁₉ - im*(Ω/2)*(a₁₅)
    da₂₀ = -im*(100*ωᵣ + Ω)*a₂₀ - im*(Ω/2)*(a₁₆)
    da₂₁ = -im*(100*ωᵣ + Ω)*a₂₁ - im*(Ω/2)*(a₁₇)
end Ω ωᵣ

# Función para extraer soluciones
function extract_solutions(solutions)
    # Obtener el número de puntos temporales y el número de variables
    n_times = length(solutions)
    n_vars = length(solutions[1])
    
    # Preallocate array para todas las soluciones
    all_sols = [Complex{Float64}[] for _ in 1:n_vars]
    
    # Extraer cada componente de la solución
    for i in 1:n_vars
        all_sols[i] = [sol[i] for sol in solutions]
    end
    
    return all_sols
end

# Función para calcular probabilidades
function calculate_probabilities(solutions)
    return [abs2.(sol) for sol in solutions]
end

# Configuración inicial
u₀ = Complex{Float64}[1; zeros(20)]  # [1, 0, 0, ..., 0]
tspan = (0.0, 0.0000005)
t = range(tspan..., length=200)

# Parámetros físicos
ω = 2π * 73.7e3  # Frecuencia para Li-6
V₀ = 50*ω

p = [V₀, ω]

# Resolver el sistema
prob = ODEProblem(my_ode10, u₀, tspan, p)
sol = solve(prob, abstol=1e-15, reltol=1e-20)

# Extraer soluciones y calcular probabilidades
solutions = extract_solutions(sol.u)
probs = calculate_probabilities(solutions)

# Plotting
function plot_probabilities(t, probs; indices=[1,4,8,12,16,20], title="Probabilidades de Transición")
    plt = plot(
        title=title,
        xlabel="Tiempo",
        ylabel="Probabilidad",
        legend=:outerright,
        dpi=1000,
        lw=2
    )
    
    labels = ["P₀(t)", "P₂(t)", "P₄(t)", "P₆(t)", "P₈(t)", "P₁₀(t)"]
    
    for (i, idx) in enumerate(indices)
        plot!(plt, t, probs[idx], label=labels[i])
    end
    
    return plt
end

# Generar gráfica
p1 = plot_probabilities(sol.t, probs, 
    title="Probabilidades para Ω = $V₀ Hz y ωᵣ = $ω Hz")

# Guardar gráfica
savefig(p1, "kapitza_dirac_simulation.png")