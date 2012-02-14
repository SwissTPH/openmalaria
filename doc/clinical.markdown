One day Case Management
=================

Sickness events
-------------------------

The case management system is entered upon sickness events, which are determined by our pathogenesis
model as a non-deterministic function of parasite density.

More specifically, the probability of a malarial event depends on parasite density. Given that such
an event occurs, there is a probability that the event is severe malaria dependent on parasite
density and a probability that a coinfection occurs dependent on age, either of which results in a
"complicated" state. Given that a malarial event does not occur, there is an age-dependent chance
that a non-malarial fever occurs.


Case management decisions
---------------------------------------------

The case management component of the model maps a sickness, age and parasite density to a treatment
type stochastically. A generalization is to view it as a decision tree with both random and
deterministic decisions. This model is only run executed when the clinical model determines an
individual to be seeking treatment, detailed elsewhere.

There are quite a few decisions that need to be made in a decision tree. Below is a list, though it
may not be complete:

    complicated/uncomplicated
    under5/over5
    
    for uncomplicated cases:
    UC1/UC2
    source
    tested
    result
    drug
    quality
    adherence
    time
    
    for complicated cases:
    source
    first test
    second test
    drug

Of these, the UC1/UC2/severe and under5/over5 selections are determined by code, as is `result`
dependant on `tested`. `first test` and `second test` only need to be known for costing. The rest
are only used in the treatment selection, so if handled generically the code doesn't need to know
about source, drug, quality, adherence and time. The generic decisions plus `tested`, `first test`
and `second test` are determined randomly based on probability tables. An example:

    <decision name="quality" depends="source,drug" values="good,bad">
        source(hospital){
            drug(SP){
                p(0.8):
                    good
                p(0.2):
                    bad
            }
            drug(AL) {
                p(0.9):
                    good
                p(0.1):
                    bad
            }
            drug(none):
                void    — we need to indicate some value to show this case wasn't forgotten (enable error checking)
        }
        source(lower level public facility){
            drug(SP){p(0.7){good}p(0.3){bad}}
            drug(AL){p(0.8){good}p(0.2){bad}}
            drug(none){void}
        }
        source(none):
            void
    </decision>


Treatment data (dosings)
----------------------------------------

Per case (complicated/uncomplicated) and per drug combination, we describe a standard dose in terms
of quantities in mg of each drug at each timepoint. Thus a standard dose for AL in an uncomplicated
case will look like:

    <treatmentData case="uncomplicated">
        <treatment drug="AL">
            <medicate drug="A" mg="80" time="0"/>
            <medicate drug="Lu" mg="480" time="0"/>
            <medicate drug="A" mg="80" time="0.333"/>
            <medicate drug="Lu" mg="480" time="0.333"/>
            <medicate drug="A" mg="80" time="1"/>
            <medicate drug="Lu" mg="480" time="1"/>
            <medicate drug="A" mg="80" time="1.5"/>
            <medicate drug="Lu" mg="480" time="1.5"/>
            <medicate drug="A" mg="80" time="2"/>
            <medicate drug="Lu" mg="480" time="2"/>
            <medicate drug="A" mg="80" time="2.5"/>
            <medicate drug="Lu" mg="480" time="2.5"/>
        </treatment>
        ....
    </treatmentData>

(where time is measured in days since the start of the case management decision making). There are
four factors which may vary a dosing schedule, for a given case and drug combination:

*   weight, actually assumed from age — multiplies quantities
*   quality — multiplies quantities for each drug, with good quality multiplying by 1.0
*   adherence — limits drug doses to a time range for each drug (for example, the range [0.3,2)
    applied to the above would remove the first and last two doses)
*   treatment seeking delay — adds a chosen number of days to each time (after adherence limits
    doses taken)

Note: we assume the patient weight from age ranges rather than making dosing accurately reflect a
patient's mass. This is to more closely reflect standard practices of prescribing a number of some
minimum dose (usually half a pill) based on age.

To deal with these variations, rather than duplicate a standard dosing schedule many times, I
propose to add the following adjustors:

*   A multiplier element for drug quantities, per drug. This would be described independently for
    each drug combination. There can be multiple multipliers with different dependencies, whos
    effects stack (e.g. a bad quality multiplier of 0.5 and age-based multiplier of 0.25 give a
    total multiplier of 0.125). An example:
    
        <treatment drug="AL">
            <medicate .../>
            <multiplier depends="quality", drugs="A,Lu">
                <case value="good">1.0,1.0</multiply>
                <case value="bad">0.5,0.5</multiply>
            </multiplier>
        </treatment>
    
    Note: sometimes this needs to depend on weight. How?
    
*   A dose-time limiter, per drug. This is also specified independently for each drug. For example:
    
        <treatment drug="AL">
            <medicate .../>
            <timeRange depends="adherence", drugs="A,Lu">
                <case value="good">0-3,0-3</case>
                <case value="miss first dose">0.3-3,0.3-3</case>
                <case value="selective">10-10,0-3</case>
            </timeRange>
        </treatment>
    
*   A time-shifter (effects applied after the dose-limiter has decided which doses to use). This is
    only intended to account for seeking delays, and so is not different for different drug
    schedules, thus is included in the base "treatmentData" element. An example, conditionally
    delaying medication of the main treatment but not of a suppository:
    
        <treatmentData case="complicated">
            ...
            <timeShift depends="source">
                <case value="prereferal and immediate referal">0</case>
                <case value="prereferal and delayed referal">1</case>
                ...
            </timeShift>
        </treatmentData>
        <treatmentData case="complicated suppositories">
            ...
            <timeShift depends="source">
                <case value="prereferal and immediate referal">0</case>
                <case value="prereferal and delayed referal">0</case>
                ...
            </timeShift>
        </treatmentData>




Severe outcomes
---------------------------

Outcomes are determined dependant on sickness events, parasite density and case management.

*   costs associated with sickness (dallies or similar plus monitary costs)
*   effects of case management on costs (hospital vs. health-centre treatment, etc.)
*   should capture parasite density effects when considering resistance


--------

Other: potential design for a probability table
---------

    <decision name="quality" depends="source,drug" values="good,bad">
        <branch depends="source" value="hospital">
            <branch depends="drug” value="SP">
                <branch depends="p" value="0.9">
                    good
                </branch>
                <branch depends="p" value="0.1">
                    bad
                </branch>
            </branch>
            <branch depends="drug” value="AL">
                <branch depends="p" value="0.8">
                    good
                </branch>
                <branch depends="p" value="0.2">
                    bad
                </branch>
            </branch>
            <branch depends="drug” value="none">
                <!-- again, for consistency, must give some value -->
                good
            </branch>
        </branch>
        <branch depends="source" value="lower level public facility">
            <branch depends="drug” value="SP">
                <branch depends="p" value="0.7">
                    good
                </branch>
                <branch depends="p" value="0.3">
                    bad
                </branch>
            </branch>
            <branch depends="drug” value="AL">
                <branch depends="p" value="0.8">
                    good
                </branch>
                <branch depends="p" value="0.2">
                    bad
                </branch>
            </branch>
            <branch depends="drug” value="none">
                <!-- again, for consistency, must give some value -->
                good
            </branch>
        </branch>
    </decision>
