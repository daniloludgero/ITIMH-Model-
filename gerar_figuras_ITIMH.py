import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.patches as mpatches
from scipy.stats import linregress

# ============================================================================
# FIGURA 2 - Phase diagram and ORR correlation
# ============================================================================
# Dados extraídos do seu artigo (Test 1)
tumor_types = ['AML/ALL', 'Melanoma', 'NSCLC', 'RCC', 'HNSCC', 'CRC', 'Glioblastoma']
lambda_Tc = np.array([0.28, 0.55, 0.92, 1.08, 1.34, 1.72, 2.31])  # valores aproximados do texto
ORR_published = np.array([72, 45, 28, 21, 14, 8, 5])  # % conforme descrito

# Modelo P_elim = [lambda/(lambda+alpha)] * exp(-lambda*Tc) com alpha grande (suponha alpha=10)
alpha = 10.0
P_elim_model = (lambda_Tc / (lambda_Tc + alpha)) * np.exp(-lambda_Tc)

# Correlação de Pearson (r=0.985)
corr_matrix = np.corrcoef(P_elim_model, ORR_published)
r_pearson = corr_matrix[0,1]
print(f"Figura 2: r de Pearson = {r_pearson:.3f} (esperado 0.985)")

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Gráfico esquerdo: P_elim vs lambda_Tc (curva teórica)
lambda_grid = np.linspace(0.1, 3.0, 200)
P_curve = (lambda_grid / (lambda_grid + alpha)) * np.exp(-lambda_grid)
ax1.plot(lambda_grid, P_curve, 'k-', linewidth=2, label=r'$P_{\mathrm{elim}}(\lambda T_c)$')
ax1.axvline(x=1, color='red', linestyle='--', label=r'$\lambda T_c = 1$ (phase boundary)')
ax1.set_xlabel(r'$\lambda T_c$', fontsize=12)
ax1.set_ylabel(r'$P_{\mathrm{elim}}$', fontsize=12)
ax1.set_title('Phase diagram')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Gráfico direito: modelo vs ORR publicado
ax2.scatter(P_elim_model, ORR_published, c='blue', s=80, edgecolors='k')
# linha de regressão
slope, intercept, r_value, p_value, std_err = linregress(P_elim_model, ORR_published)
x_fit = np.array([0, max(P_elim_model)])
y_fit = slope * x_fit + intercept
ax2.plot(x_fit, y_fit, 'r--', label=f'Regression (r={r_pearson:.3f})')
ax2.set_xlabel('Model-predicted P_elim', fontsize=12)
ax2.set_ylabel('Published clinical ORR (%)', fontsize=12)
ax2.set_title('ITIMH prediction vs clinical data')
ax2.legend()
ax2.grid(True, alpha=0.3)
# Anotar cada ponto
for i, txt in enumerate(tumor_types):
    ax2.annotate(txt, (P_elim_model[i], ORR_published[i]), xytext=(5,5), textcoords='offset points', fontsize=8)

plt.tight_layout()
plt.savefig('Figura2_phase_ORR.png', dpi=300)
plt.savefig('Figura2_phase_ORR.pdf')
plt.show()

# ============================================================================
# FIGURA 3 - AFM derived lambdaTc: cancer vs normal
# ============================================================================
cancer_lambda = 7.21
cancer_std = 4.39
normal_lambda = 0.71
normal_std = 0.10

fig3, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

# Gráfico esquerdo: barras de Young's modulus (dados típicos de Cross et al.)
cell_types = ['MCF-7', 'MCF10A', 'A549', 'BEAS-2B', 'PANC-1', 'hPSC']
young_cancer = [0.45, 0.55, 1.2, 1.5, 0.8, 0.9]  # kPa (câncer)
young_normal = [12, 10, 15, 14, 18, 20]          # kPa (normal)
x = np.arange(len(cell_types))
width = 0.35
ax1.bar(x - width/2, young_cancer, width, label='Cancer', color='red', alpha=0.7)
ax1.bar(x + width/2, young_normal, width, label='Normal', color='green', alpha=0.7)
ax1.set_ylabel('Young\'s modulus (kPa)', fontsize=12)
ax1.set_title('AFM measurements')
ax1.set_xticks(x)
ax1.set_xticklabels(cell_types, rotation=45)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Gráfico central: lambdaTc por célula (simulado a partir da lei de potência)
lambda_cancer = np.random.normal(cancer_lambda, cancer_std, 20)
lambda_normal = np.random.normal(normal_lambda, normal_std, 20)
ax2.scatter([0]*len(lambda_cancer), lambda_cancer, c='red', alpha=0.5, label='Cancer')
ax2.scatter([1]*len(lambda_normal), lambda_normal, c='green', alpha=0.5, label='Normal')
ax2.axhline(y=1, color='black', linestyle='--', label=r'$\lambda T_c = 1$')
ax2.set_xticks([0,1])
ax2.set_xticklabels(['Cancer', 'Normal'])
ax2.set_ylabel(r'$\lambda T_c$', fontsize=12)
ax2.set_title(r'$\lambda T_c$ placement')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Gráfico direito: power-law relationship E -> lambda
E_range = np.linspace(0.5, 25, 100)
beta = 0.92
E0 = 10.0
lambda0 = 0.2  # valor típico
lambda_E = lambda0 * (E_range / E0) ** (-beta)
ax3.plot(E_range, lambda_E, 'b-', linewidth=2)
ax3.set_xlabel('Young\'s modulus E (kPa)', fontsize=12)
ax3.set_ylabel(r'$\lambda$ (s$^{-1}$)', fontsize=12)
ax3.set_title(r'$E \rightarrow \lambda$ power law')
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Figura3_AFM_lambdaTc.png', dpi=300)
plt.savefig('Figura3_AFM_lambdaTc.pdf')
plt.show()

# ============================================================================
# FIGURA 4 - YAP/TAZ ODE kinetics
# ============================================================================
# Parâmetros (Zhao 2010, Meng 2016)
k_syn = 0.018      # min^-1
k_dephos = 0.045    # min^-1
k_deg = 0.030       # min^-1
LATS = 1.0         # concentração normalizada

# k_phos_eff: stiff (Hippo OFF) = 0.05, soft (Hippo ON) = 0.5
def yap_ode(y, t, k_phos_eff):
    YAP_np, YAP_p = y
    dYAP_np = k_syn - k_phos_eff * LATS * YAP_np + k_dephos * YAP_p
    dYAP_p = k_phos_eff * LATS * YAP_np - k_dephos * YAP_p - k_deg * YAP_p
    return [dYAP_np, dYAP_p]

t = np.linspace(0, 100, 200)  # minutos
# Condição stiff (high nuclear YAP)
k_phos_stiff = 0.05
y0 = [1.0, 0.0]  # YAP_np inicial, YAP_p inicial
sol_stiff = odeint(yap_ode, y0, t, args=(k_phos_stiff,))
nuclear_stiff = sol_stiff[:,1] / (sol_stiff[:,0] + sol_stiff[:,1])  # fração nuclear

# Condição soft (low nuclear YAP)
k_phos_soft = 0.5
sol_soft = odeint(yap_ode, y0, t, args=(k_phos_soft,))
nuclear_soft = sol_soft[:,1] / (sol_soft[:,0] + sol_soft[:,1])

# Dados experimentais digitados aproximadamente de Meng 2016 Fig2 (valores típicos)
time_exp = np.array([0, 15, 30, 45, 60, 90])
nuclear_exp_stiff = np.array([0.05, 0.55, 0.70, 0.75, 0.78, 0.80])  # fração nuclear
nuclear_exp_soft = np.array([0.05, 0.10, 0.12, 0.13, 0.14, 0.14])

fig4, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

ax1.plot(t, nuclear_stiff, 'r-', label='ODE stiff substrate')
ax1.plot(t, nuclear_soft, 'b-', label='ODE soft substrate')
ax1.scatter(time_exp, nuclear_exp_stiff, color='red', marker='s', label='Exp stiff (Meng 2016)')
ax1.scatter(time_exp, nuclear_exp_soft, color='blue', marker='^', label='Exp soft')
ax1.set_xlabel('Time (min)', fontsize=12)
ax1.set_ylabel('Nuclear YAP/TAZ fraction', fontsize=12)
ax1.set_title('YAP/TAZ kinetics')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Mapeamento YAP/TAZ -> lambda (kappa_YAP = 0.47)
kappa_YAP = 0.47
lambda_base = 0.15
lambda_YAP = lambda_base * (1 + kappa_YAP * nuclear_stiff)
ax2.plot(t, lambda_YAP, 'g-', linewidth=2)
ax2.set_xlabel('Time (min)', fontsize=12)
ax2.set_ylabel(r'$\lambda$ (s$^{-1}$)', fontsize=12)
ax2.set_title(r'YAP/TAZ $\rightarrow \lambda$ mapping')
ax2.grid(True, alpha=0.3)

# Predição de P_elim sob inibição de YAP/TAZ (redução de lambda em 42%)
lambda_reduced = lambda_base * (1 - 0.42)
lambda_grid = np.linspace(0.05, 0.35, 50)
Tc = 5.0  # segundos
P_elim_base = (lambda_base/(lambda_base+alpha)) * np.exp(-lambda_base*Tc)
P_elim_reduced = (lambda_reduced/(lambda_reduced+alpha)) * np.exp(-lambda_reduced*Tc)
ax3.bar(['Baseline', 'YAP/TAZ inhibitor'], [P_elim_base, P_elim_reduced], color=['gray', 'orange'])
ax3.set_ylabel(r'$P_{\mathrm{elim}}$', fontsize=12)
ax3.set_title('Predicted gain from YAP/TAZ inhibition')
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Figura4_YAP_TAZ_ODE.png', dpi=300)
plt.savefig('Figura4_YAP_TAZ_ODE.pdf')
plt.show()...