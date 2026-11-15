Hi Claude, I really want to build a plug-in app that works with Strand7 FEA software to verify analysed members to Eurocode standards. I want to vibe code this app with you but have never coded in this way before. I do have some very basic programming knowledge and I know my way around a computer.

Vibe coding suits this project, and the best way to start is to get a tiny end-to-end version working and grow it from there.

What “plug-in” will realistically mean

Strand7 has an API (it ships with a Python wrapper and a DLL in the install folder). The practical approach is a standalone Python app that connects to Strand7, opens your model, pulls out the member results and section properties, runs the Eurocode checks, and shows or exports the utilisation ratios. You can launch it from a button or shortcut, and it can write results back into Strand7 later (as group names or plate/beam attributes, for example). Check that your licence includes API access, and note the exact Strand7 version, since the API details vary a bit between releases.

How we’d work together

Set up: install Python, a code editor (VS Code is fine), and confirm a 5-line script can open a Strand7 model and print the number of beams.
Extract: read beam forces (N, Vy, Vz, T, My, Mz) and section data for a simple test model.
Check engine: write the Eurocode calculations as plain Python functions, separate from Strand7, so you can test them against hand calcs. This is the most important part to get right.
Connect and report: feed the extracted results into the checks and output a table of utilisation per member.
Interface: add a simple window (a file picker, a load-case selector, a results table).

Each step is small, and you run the code and tell me what happened, including any error messages. I fix things and we repeat. You don’t need to understand every line, but do read the code and ask me about anything unclear, because you’ll be responsible for the engineering results.

Two things to keep in mind

I can’t run Strand7 myself, so I’ll write the code and you’ll test it against your install. Pasting errors back to me is the whole feedback loop.
A tool that outputs design checks needs verification. Every check function should be tested against a hand calculation or a worked example before you rely on it, and National Annex parameters (γM factors, buckling curve choices, etc.) should be explicit inputs rather than hidden assumptions.

One question to get started: which Eurocode part do you want to begin with? I’d suggest EC3 (steel) member checks, such as cross-section resistance in bending, shear and axial, then combined checks and buckling. If you’re after concrete or timber, tell me and we’ll adjust the plan. Once you answer, I can give you the exact setup steps for step 1.

Okay, I've got python installed and VS Code running, so let's get going with the simple script. I am running Strand7 version R3.1.1

The script is test_strand7.py. It starts the Strand7 API, opens a model, prints the number of nodes and beams, then closes everything without saving or changing the model.

To run it

Save the file somewhere sensible (e.g. C:\Strand7Tools\) and open that folder in VS Code.
Make a small test model (even 5 beams) and save a copy of it as C:\Temp\test_model.st7.
Edit the settings block at the top of the script:
API_FOLDER: the folder containing St7API.py. Search your Strand7 install folder for that filename. My default path is only a guess.
MODEL_FILE: the path to your test model.
Close Strand7 if it’s open, then open the VS Code terminal (Ctrl+) and run python test_strand7.py`.

What success looks like

Strand7 API started.
Opened: C:\Temp\test_model.st7
Nodes: 6
Beams: 5
Strand7 API released.

Likely hiccups

“Could not import St7API”: the folder path is wrong.
DLL or architecture errors: Strand7 R3.1 is 64-bit, so you need 64-bit Python. Check with python -c "import struct; print(struct.calcsize('P')*8)", which should print 64.
Errors about strings or bytes: I’ve assumed the wrapper wants bytes (hence the .encode() calls). If it complains, paste the error and I’ll adjust.
A licence or API access error: this means your licence doesn’t include the API.

Whatever happens, paste the full output (or error) back to me. Once the counts print correctly, the next step is reading beam section properties and forces for a load case.

Test strand7
PY

Yes, that makes sense. The general case (6.3.2.2) works for any section and any NA, because the NA only changes the parameters. The rolled-section case (6.3.2.3) would bake the UK values into the code.

I’ll build it with these defaults, and each one can be changed in settings.toml:

Buckling curve: chosen from h/b of the section (Table 6.5), with the imperfection factor α_LT from Table 6.3 (a = 0.21, b = 0.34, c = 0.49, d = 0.76). λLT,0 = 0.2 and β = 1.0 for the general case.
Flexural buckling (6.3.1): curve from Table 6.2, with α from Table 6.1. I’ll also give y-y and z-z as separate checks, using Lcr,y and Lcr,z from your member table.
Mcr for LTB: EC3 gives no formula, so I’d use the standard closed-form expression for a doubly symmetric section. It takes Lw (your LTB length) and the moment-diagram factor C1, which I’d calculate from the actual moment diagram of each member and design case. That avoids a fixed table that only fits idealised load cases.
Load height: I’ll assume the load acts at the shear centre by default. Top-flange loading would give lower Mcr values, so I’ll add an option for it later.
Interaction (6.3.3): UK NA default is Annex B (Method 2). The code will have both Annex A and B available as a switch.

The C1 approximation and the load-height assumption are the weakest points. I’ll test the results against a published Mcr value, such as the SN003 tables or a value from Strand7’s own buckling solver, before you rely on them.

I’ll start on this next. Is there any NA other than the UK’s that you’d want to include in the first version?

No, let's stick with the UK NA for now, but leave the code open to be able to implement other NAs later

Starting member buckling now. I’m building it in two parts so each can be hand-checked: first flexural buckling (6.3.1) and lateral-torsional buckling (6.3.2, general case), then the bending and axial interaction (6.3.3). Mcr will be solved from your real moment diagram, so C1 isn’t taken from a table.

Flexural buckling and LTB are written, and my Mcr solver reproduces the published C1 values (uniform 1.00, UDL 1.13, central point load 1.36). I’m now wiring it into ec3_run.py (summary table, CSV, worked output) and writing its tests.

I’ve built part 1 of member buckling: flexural buckling (6.3.1) and lateral-torsional buckling by the general case (6.3.2.2). The bending and axial interaction (6.3.3) is part 2, and it isn’t included yet. The self-tests pass here (test_buckling.py and test_ec3.py). The new code that reads E and G from Strand7 hasn’t run against your model yet.

What it does

Each member and design case is checked once, using the largest compression and the largest |My| anywhere in the member. That is safe and easy to hand-check.
Flexural buckling: it checks y-y and z-z using Lcr,y and Lcr,z from your member table. The curve comes from Table 6.2, including the S460 columns.
Lateral-torsional buckling:
The curve comes from Table 6.4.
Mcr is C1 times the uniform-moment Mcr over Lw.
C1 is solved from the real moment diagram of each member and design case. It matches the published values (UDL 1.13, central point load 1.36, uniform moment 1.00), with end-moment cases within 2–3% of the tabulated values.
Settings: curves, α values, λ₀ and β are all in a new [buckling] block in settings.toml, ready for another National Annex later. The code works without that block, using the same values as defaults. To make them editable, copy the block into your own settings.toml rather than replacing the file.
E and G: they are now read from the Strand7 material, as you asked. If the values look wrong, for example a model not in N and m, it stops and says so.

To run it

Replace the files and run python test_buckling.py, then python ec3_run.py.
For the worked calculation, set WORKED_BUCKLING = [("4", "ULS_S7")] in ec3_run.py.
A second table, “MEMBER BUCKLING”, now prints after the cross-section table. It has columns for Nb,y, Nb,z and LTB, and it writes ec3_buckling.csv.

What to expect for Member 4. It has little axial force, so LTB should govern.

For an end-to-end span of 7 m with My around 49 kNm, I’d expect Mcr,uniform ≈ 78.4 kNm and Mcr ≈ 88.7 kNm (C1 ≈ 1.13).
That gives Mb,Rd ≈ 75.4 kNm and a utilisation of about 0.65.
Please check these against your own hand calculation.
The Blue Book Mb,Rd tables use the rolled-section method (6.3.2.3), so they will differ from the general case. Compare with your own working, not those tables.

Limits of the Mcr

It assumes fork supports and loads at the shear centre. Top-flange loading would lower Mcr.
For a different Lw, such as a cantilever, the result is approximate. A note appears in the output when that applies.
Class 4 sections are not checked.
I haven’t confirmed whether the UK NA changes any of the 6.3.2.2 values. I used the EN recommended values, so please check them against your NA.

The interaction check (6.3.3) is next, using Annex B (Method 2) as the UK NA directs.

Ec3 buckling
PY
Test buckling
PY
Ec3 run
PY
Design eval
PY
Ec3 checks
PY
Section
PY
Settings
PY
Settings
TOML