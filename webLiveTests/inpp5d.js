import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://trena.systemsbiology.net`;

test('launch INPP5D app', async function foo(t){

    const launchAppLink = Selector("a").withText("trena INPP5D");
    const launchingNotification = Selector("span").withText("trena INPP5D");
    const igvPresentInIframe = Selector("body > div > h2").withText("INPP5D Trena Model & Disruptions");
    const igvCurrentGenomeDiv = Selector("div.igv-current_genome");
    //const igvCurrentGenomeDiv = Selector("#igv-current_genome");

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
