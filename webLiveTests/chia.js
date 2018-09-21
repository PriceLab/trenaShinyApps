import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://chia.systemsbiology.net`;

test('launch chia app', async function foo(t){

    const pageLoaded = Selector("#selectDatasetMenuLabel");
    await t
        .expect(pageLoaded.exists).ok()
    });
