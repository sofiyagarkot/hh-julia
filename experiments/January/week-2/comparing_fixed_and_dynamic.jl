include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = true

evaltspan = (0.0, 50.0)
dts = 10.0 .^ range(-2, -3, length=10)[begin:end-1]
DM = FixedDiffusion()
algorithms = [
    EK0(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
    EK0(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
#      EK1(prior=IWP(5), smooth=SMOOTH, diffusionmodel=DM)
    ]

names = [  
"IWP(1), FixedDiffusion", 
"IWP(2), FixedDiffusion", 
"IWP(3), FixedDiffusion", 
"IWP(4), FixedDiffusion", 
# "IWP(5)"
]

algorithms2 = [
    EK0(prior=IWP(1), smooth=SMOOTH),
    EK0(prior=IWP(2), smooth=SMOOTH),
    EK0(prior=IWP(3), smooth=SMOOTH), 
    EK0(prior=IWP(4), smooth=SMOOTH), 
    EK1(prior=IWP(5), smooth=SMOOTH)
    ]

names2 = [  
    "IWP(1), DynamicDiffusion", 
    "IWP(2), DynamicDiffusion", 
    "IWP(3), DynamicDiffusion", 
    "IWP(4), DynamicDiffusion", 
    # "IWP(5)"
    ]


prob = define_problem()


plot_2_work_precision_plots(
    names::Array, 
    algorithms :: Array, 
    names2::Array, 
    algorithms2 :: Array, 
    prob::ODEProblem; 

    DENSE = DENSE, 
    SAVE_EVERYSTEP = false, 
    to_save=TO_SAVE, 
    to_save_path="./visuals/January/diffusions/fixed_and_dynamic_wp_EK0_IWP.png",
    title="EK0, Fixed diffusion vs. Dynamic diffusion with adaptive",
    adaptive=ADAPTIVE,
    dts=dts
    )

# EK1 
algorithms = [
    EK1(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
    EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
    EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
    EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
#      EK1(prior=IWP(5), smooth=SMOOTH, diffusionmodel=DM)
    ]

names = [  
"IWP(1), FixedDiffusion", 
"IWP(2), FixedDiffusion", 
"IWP(3), FixedDiffusion", 
"IWP(4), FixedDiffusion", 
# "IWP(5)"
]

algorithms2 = [
    EK1(prior=IWP(1), smooth=SMOOTH),
    EK1(prior=IWP(2), smooth=SMOOTH),
    EK1(prior=IWP(3), smooth=SMOOTH), 
    EK1(prior=IWP(4), smooth=SMOOTH), 
    EK1(prior=IWP(5), smooth=SMOOTH)
    ]

names2 = [  
    "IWP(1), DynamicDiffusion", 
    "IWP(2), DynamicDiffusion", 
    "IWP(3), DynamicDiffusion", 
    "IWP(4), DynamicDiffusion", 
    "IWP(5), DynamicDiffusion"
    ]


prob = define_problem()


plot_2_work_precision_plots(
    names::Array, 
    algorithms :: Array, 
    names2::Array, 
    algorithms2 :: Array, 
    prob::ODEProblem; 

    DENSE = DENSE, 
    SAVE_EVERYSTEP = false, 
    to_save=TO_SAVE, 
    to_save_path="./visuals/January/diffusions/fixed_and_dynamic_wp_EK1_IWP.png",
    title="EK1, Fixed diffusion vs. Dynamic diffusion with adaptive",
    adaptive=ADAPTIVE,
    dts=dts
    )