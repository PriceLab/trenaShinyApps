import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://trena.systemsbiology.net`;

test('launch FRD3 app', async function foo(t){

    const launchAppLink = Selector("a").withText("trena FRD3");
    const launchingNotification = Selector("span").withText("trena FRD3");
    const igvPresentInIframe = Selector("h2").withText("FRD3 Trena Models");
    const igvCurrentGenomeDiv = Selector("#igv-current_genome");

    await t
        .typeText('#username', 'paul')
        .typeText('#password', 'leo.paul')
        .click('.form-signin > button')
        .expect(launchAppLink.exists).ok()
        .click(launchAppLink)
        .wait(5000)
        .switchToIframe("#shinyframe")
        .expect(igvPresentInIframe.exists).ok()
        .expect(igvCurrentGenomeDiv.exists).ok();
    });
