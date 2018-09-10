import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://trena.systemsbiology.net`;

test('launch IREG1 app', async function foo(t){

    const launchAppLink = Selector("a").withText("trena IREG1");
    const launchingNotification = Selector("span").withText("trena IREG1");
    const igvPresentInIframe = Selector("body > div > h2").withText("IREG1 Trena Models");
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
